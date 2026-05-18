#!/usr/bin/env python3
"""HTTP server for the sediment-sim bridge.
Accepts sim state as JSON, runs OpenFOAM, returns velocity field.

POST /bake  — {"obstacles": [...], "inflow": 2.5} → runs solver → {"ux": [...], "uy": [...]}
GET /status — returns current job progress
"""

import json
import subprocess
import sys
import tempfile
import threading
import time
import re
from pathlib import Path
from http.server import HTTPServer, BaseHTTPRequestHandler

# Add parent to path for case_generator import
sys.path.insert(0, str(Path(__file__).parent))
from case_generator import SimState, Obstacle, create_case, generate_initial_U

PORT = 8090
CASE_DIR = Path(__file__).parent / "bake_case"


class JobManager:
    """Manages background OpenFOAM jobs."""
    def __init__(self):
        self.current_job = None
        self.lock = threading.Lock()

    def start_bake(self, state: dict) -> str:
        """Start a bake job from sim state dict. Returns job ID."""
        with self.lock:
            job_id = f"bake_{int(time.time())}"
            self.current_job = {
                "id": job_id,
                "status": "starting",
                "progress": 0,
                "result": None,
                "error": None,
                "started": time.time()
            }
        
        thread = threading.Thread(target=self._run_bake, args=(job_id, state))
        thread.daemon = True
        thread.start()
        return job_id

    def _run_bake(self, job_id: str, state_dict: dict):
        """Run the full OpenFOAM pipeline in background."""
        try:
            self._update(job_id, "meshing", 10)
            
            # Build sim state from JSON
            obstacles = []
            for ob in state_dict.get("obstacles", []):
                obstacles.append(Obstacle(
                    type=ob["type"],
                    x=ob["x"],
                    y=ob["y"],
                    size=ob.get("size", 40),
                    rotation=ob.get("rotation", 0)
                ))
            state = SimState(
                obstacles=obstacles,
                inflow_velocity=state_dict.get("inflow", 2.5),
                width_px=state_dict.get("width", 1280),
                height_px=state_dict.get("height", 577)
            )
            
            # Clean and regenerate case
            if CASE_DIR.exists():
                import shutil
                shutil.rmtree(CASE_DIR)
            create_case(state, CASE_DIR)
            
            of_prefix = "source /usr/share/openfoam/etc/bashrc && "
            
            # blockMesh
            self._update(job_id, "blockMesh", 20)
            proc = subprocess.run(
                of_prefix + f"blockMesh -case {CASE_DIR}",
                shell=True, executable='/bin/bash',
                capture_output=True, text=True, timeout=300
            )
            if proc.returncode != 0:
                raise RuntimeError(f"blockMesh failed: {proc.stderr[-300:]}")
            
            # snappyHexMesh
            self._update(job_id, "snappyHexMesh", 40)
            proc = subprocess.run(
                of_prefix + f"snappyHexMesh -overwrite -case {CASE_DIR}",
                shell=True, executable='/bin/bash',
                capture_output=True, text=True, timeout=600
            )
            if proc.returncode != 0:
                raise RuntimeError(f"snappyHexMesh failed: {proc.stderr[-300:]}")
            
            # Rebuild 0/ fields with obstacle patches (after snappy creates them)
            (CASE_DIR / "0" / "U").write_text(generate_initial_U(state))
            from case_generator import generate_initial_p, generate_k, generate_epsilon, generate_nut
            (CASE_DIR / "0" / "p").write_text(generate_initial_p(state))
            (CASE_DIR / "0" / "k").write_text(generate_k(state))
            (CASE_DIR / "0" / "epsilon").write_text(generate_epsilon(state))
            (CASE_DIR / "0" / "nut").write_text(generate_nut(state))
            
            # simpleFoam with progress parsing
            self._update(job_id, "simpleFoam", 50)
            proc = subprocess.Popen(
                of_prefix + f"simpleFoam -case {CASE_DIR}",
                shell=True, executable='/bin/bash',
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            
            last_time = 0
            for line in proc.stdout:
                # Parse "Time = N" from output
                m = re.search(r"Time\s*=\s*(\d+)", line)
                if m:
                    t = int(m.group(1))
                    if t != last_time:
                        last_time = t
                        # Progress from 50 to 95, scaled by iteration count
                        pct = min(50 + int(45 * t / 500), 95)
                        self._update(job_id, f"iter {t}", pct)
            
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"simpleFoam failed with exit code {proc.returncode}")
            
            # Parse velocity field
            self._update(job_id, "parsing results", 98)
            velocity = self._parse_velocity()
            
            self._update(job_id, "complete", 100, result=velocity)
            
        except Exception as e:
            self._update(job_id, "failed", 100, error=str(e))

    def _parse_velocity(self):
        """Parse U field and return simplified velocity data."""
        # Find latest time directory
        time_dirs = sorted(
            [d for d in CASE_DIR.iterdir() 
             if d.is_dir() and d.name.replace('.', '').replace('-', '').isdigit()],
            key=lambda d: float(d.name)
        )
        if not time_dirs:
            return None
        
        u_file = time_dirs[-1] / "U"
        if not u_file.exists():
            return None
        
        content = u_file.read_text()
        start = content.find("internalField")
        section = content[start:].split('\n')
        
        vals = []
        count = None
        reading = False
        
        for line in section:
            s = line.strip()
            if not reading:
                if 'nonuniform' in s:
                    continue
                if count is None and s.isdigit():
                    count = int(s)
                    continue
                if count is not None and '(' in s:
                    reading = True
                    # Parse first value if on same line
                    idx = s.find('(')
                    rest = s[idx+1:]
                    if rest and ')' not in rest:
                        nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', rest)
                        if len(nums) >= 3:
                            vals.append((float(nums[0]), float(nums[1])))
                    continue
            else:
                if s == ')' or s.startswith(');'):
                    break
                if s and not s.startswith('//'):
                    nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', s)
                    if len(nums) >= 3:
                        vals.append((float(nums[0]), float(nums[1])))
        
        # Downsample for transfer: take every Nth value to get ~2000 points
        step = max(1, len(vals) // 2000)
        sampled = vals[::step]
        
        return {
            "n_cells": len(vals),
            "sampled": len(sampled),
            "ux": [v[0] for v in sampled],
            "uy": [v[1] for v in sampled],
            "ux_min": min(v[0] for v in vals),
            "ux_max": max(v[0] for v in vals),
            "ux_mean": sum(v[0] for v in vals) / len(vals),
            "reverse_cells": sum(1 for v in vals if v[0] < 0),
        }

    def _update(self, job_id, status, progress, result=None, error=None):
        with self.lock:
            if self.current_job and self.current_job["id"] == job_id:
                self.current_job["status"] = status
                self.current_job["progress"] = progress
                if result is not None:
                    self.current_job["result"] = result
                if error is not None:
                    self.current_job["error"] = error

    def get_status(self):
        with self.lock:
            if self.current_job is None:
                return {"status": "idle"}
            j = self.current_job.copy()
            j["elapsed"] = time.time() - j["started"]
            return j


job_manager = JobManager()


class BridgeHandler(BaseHTTPRequestHandler):
    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        if self.path == "/status":
            status = job_manager.get_status()
            self._send_json(status)
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        if self.path == "/bake":
            content_length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_length)
            try:
                state = json.loads(body)
            except json.JSONDecodeError:
                self.send_response(400)
                self._send_json({"error": "Invalid JSON"})
                return
            
            job_id = job_manager.start_bake(state)
            self._send_json({"job_id": job_id, "status": "started"})
        else:
            self.send_response(404)
            self.end_headers()

    def _send_json(self, data):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(data).encode())

    def log_message(self, format, *args):
        print(f"[{time.strftime('%H:%M:%S')}] {args[0]}")


def main():
    print(f"Bridge server starting on http://localhost:{PORT}")
    print("Endpoints: POST /bake  |  GET /status")
    server = HTTPServer(("0.0.0.0", PORT), BridgeHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()
