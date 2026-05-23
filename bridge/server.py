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
            
            # Check for bed morphology
            depo_png = state_dict.get("depo_png", "")
            bed_profile = state_dict.get("bed_profile", [])
            has_bed = len(bed_profile) >= 2  # Need at least 2 points for a profile
            
            of_prefix = "source /usr/share/openfoam/etc/bashrc && "
            
            if has_bed:
                # Standard blockMesh first, then displace bottom vertices
                self._update(job_id, "blockMesh", 15)
                proc = subprocess.run(
                    of_prefix + f"blockMesh -case {CASE_DIR}",
                    shell=True, executable='/bin/bash',
                    capture_output=True, text=True, timeout=300
                )
                if proc.returncode != 0:
                    raise RuntimeError(f"blockMesh failed: {proc.stderr[-300:]}")
                
                # Displace bottom by bed heightfield
                self._update(job_id, "displacing bed", 20)
                from displace_bed import displace_bottom
                displace_bottom(CASE_DIR, depo_png, bed_profile=bed_profile)
            else:
                # Standard blockMesh with flat bottom
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
            
            # Note: field files from create_case already have obstacle patches
            # that match the SimState. If snappy creates extra patches, 
            # they'll need to be handled. For now this works for our obstacle set.
            
            # simpleFoam with progress parsing
            self._update(job_id, "simpleFoam", 50)
            proc = subprocess.Popen(
                of_prefix + f"simpleFoam -case {CASE_DIR}",
                shell=True, executable='/bin/bash',
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            
            last_time = 0
            try:
                for line in proc.stdout:
                    m = re.search(r"Time\s*=\s*(\d+)", line)
                    if m:
                        t = int(m.group(1))
                        if t != last_time:
                            last_time = t
                            pct = min(50 + int(45 * t / 500), 95)
                            self._update(job_id, f"iter {t}", pct)
            except Exception:
                pass
            
            # Ensure process finishes and check exit code
            proc.wait(timeout=600)
            if proc.returncode != 0:
                raise RuntimeError(f"simpleFoam failed with exit code {proc.returncode}")
            
            # Parse velocity field
            self._update(job_id, "parsing results", 98)
            velocity = self._parse_velocity()
            
            self._update(job_id, "complete", 100, result=velocity)
            
        except Exception as e:
            self._update(job_id, "failed", 100, error=str(e))

    def _parse_velocity(self):
        """Parse U field and return velocity data on a grid matching Stam solver."""
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
        
        # Interpolate to Stam grid (160×90) using actual cell centers
        from scipy.spatial import cKDTree
        from cell_centers import compute_cell_centers
        from case_generator import DOMAIN_W, DOMAIN_H
        STAM_NX, STAM_NY = 160, 90
        
        # Get actual cell center positions from the mesh
        cell_centers = compute_cell_centers(CASE_DIR / "constant" / "polyMesh")
        
        # Build KD-tree from cell centers (normalized to [0,1])
        positions = []
        for cx, cy in cell_centers:
            positions.append((cx / DOMAIN_W, cy / DOMAIN_H))
        
        # If cell count mismatch (e.g., snappy changed cell count), 
        # pad or truncate velocity values
        n_vel = len(vals)
        n_cells = len(positions)
        if n_vel != n_cells:
            # Truncate to smaller count
            n = min(n_vel, n_cells)
            vals = vals[:n]
            positions = positions[:n]
        
        tree = cKDTree(positions)
        
        # Query Stam grid points
        mag_flat = []
        for j in range(STAM_NY):
            for i in range(STAM_NX):
                qx = (i + 0.5) / STAM_NX
                qy = (STAM_NY - 1 - j + 0.5) / STAM_NY  # Flip Y
                dist, idx = tree.query((qx, qy), k=1)
                ux, uy = vals[idx]
                mag_flat.append((ux*ux + uy*uy) ** 0.5)
        
        return {
            "n_cells": len(vals),
            "ux_min": min(v[0] for v in vals),
            "ux_max": max(v[0] for v in vals),
            "ux_mean": sum(v[0] for v in vals) / len(vals),
            "reverse_cells": sum(1 for v in vals if v[0] < 0),
            # Grid data for side-by-side comparison (160×90 = Stam resolution)
            "grid": {
                "nx": STAM_NX, "ny": STAM_NY,
                "mag": mag_flat,  # Flattened row-major
                "mag_max": max(mag_flat) if mag_flat else 1.0,
            }
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


# ── Post-snappy field generators (handles arbitrary patch names) ──

def _parse_boundary_patches(boundary_file):
    """Read patch names from a polyMesh/boundary file."""
    text = boundary_file.read_text()
    patches = []
    for m in re.finditer(r'^\s+(\w+)\s*$', text, re.MULTILINE):
        name = m.group(1)
        # Filter out OpenFOAM metadata keywords, but keep actual patch names
        if name not in ('{', '}', 'type', 'inGroups', 'nFaces', 'startFace', 
                         'patch', 'wall', 'empty', 'symmetry', 'cyclic',
                         'nonuniform', 'uniform', 'calculated', 'fixedValue',
                         'zeroGradient', 'slip', 'kqRWallFunction', 
                         'epsilonWallFunction', 'nutkWallFunction',
                         'class', 'format', 'version', 'object', 'location'):
            if name not in patches and not name.startswith('//') and not name.startswith('#'):
                patches.append(name)
    return patches


def _generate_U_with_patches(state, patch_names):
    u_in = state.inflow_velocity
    bc_entries = []
    for p in patch_names:
        if p == 'inlet':
            bc_entries.append(f"    inlet\n    {{\n        type            fixedValue;\n        value           uniform ({u_in} 0 0);\n    }}")
        elif p == 'outlet':
            bc_entries.append(f"    outlet\n    {{\n        type            zeroGradient;\n    }}")
        elif p in ('top', 'bottom'):
            bc_entries.append(f"    {p}\n    {{\n        type            slip;\n    }}")
        elif p == 'frontAndBack':
            bc_entries.append(f"    frontAndBack\n    {{\n        type            empty;\n    }}")
        else:
            bc_entries.append(f"    {p}\n    {{\n        type            fixedValue;\n        value           uniform (0 0 0);\n    }}")
    
    bc_block = '\n'.join(bc_entries)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({u_in} 0 0);

boundaryField
{{
{bc_block}
}}
"""


def _generate_p_with_patches(state, patch_names):
    bc_entries = []
    for p in patch_names:
        if p == 'inlet':
            bc_entries.append(f"    inlet\n    {{\n        type            zeroGradient;\n    }}")
        elif p == 'outlet':
            bc_entries.append(f"    outlet\n    {{\n        type            fixedValue;\n        value           uniform 0;\n    }}")
        elif p in ('top', 'bottom'):
            bc_entries.append(f"    {p}\n    {{\n        type            slip;\n    }}")
        elif p == 'frontAndBack':
            bc_entries.append(f"    frontAndBack\n    {{\n        type            empty;\n    }}")
        else:
            bc_entries.append(f"    {p}\n    {{\n        type            zeroGradient;\n    }}")
    bc_block = '\n'.join(bc_entries)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
{bc_block}
}}
"""


def _generate_k_with_patches(state, patch_names):
    bc_entries = []
    for p in patch_names:
        if p == 'inlet':
            bc_entries.append(f"    inlet   {{ type fixedValue; value uniform 0.001; }}")
        elif p == 'outlet':
            bc_entries.append(f"    outlet  {{ type zeroGradient; }}")
        elif p in ('top', 'bottom'):
            bc_entries.append(f"    {p}     {{ type kqRWallFunction; value uniform 0.001; }}")
        elif p == 'frontAndBack':
            bc_entries.append(f"    frontAndBack {{ type empty; }}")
        else:
            bc_entries.append(f"    {p}     {{ type kqRWallFunction; value uniform 0.001; }}")
    bc_block = '\n'.join(bc_entries)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.001;

boundaryField
{{
{bc_block}
}}
"""


def _generate_eps_with_patches(state, patch_names):
    bc_entries = []
    for p in patch_names:
        if p == 'inlet':
            bc_entries.append(f"    inlet   {{ type fixedValue; value uniform 0.0001; }}")
        elif p == 'outlet':
            bc_entries.append(f"    outlet  {{ type zeroGradient; }}")
        elif p in ('top', 'bottom'):
            bc_entries.append(f"    {p}     {{ type epsilonWallFunction; value uniform 0.0001; }}")
        elif p == 'frontAndBack':
            bc_entries.append(f"    frontAndBack {{ type empty; }}")
        else:
            bc_entries.append(f"    {p}     {{ type epsilonWallFunction; value uniform 0.0001; }}")
    bc_block = '\n'.join(bc_entries)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 0.0001;

boundaryField
{{
{bc_block}
}}
"""


def _generate_nut_with_patches(state, patch_names):
    bc_entries = []
    for p in patch_names:
        if p == 'inlet':
            bc_entries.append(f"    inlet   {{ type calculated; value uniform 0; }}")
        elif p == 'outlet':
            bc_entries.append(f"    outlet  {{ type calculated; value uniform 0; }}")
        elif p in ('top', 'bottom'):
            bc_entries.append(f"    {p}     {{ type nutkWallFunction; value uniform 0; }}")
        elif p == 'frontAndBack':
            bc_entries.append(f"    frontAndBack {{ type empty; }}")
        else:
            bc_entries.append(f"    {p}     {{ type nutkWallFunction; value uniform 0; }}")
    bc_block = '\n'.join(bc_entries)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
{bc_block}
}}
"""


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
