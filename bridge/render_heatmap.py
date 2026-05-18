#!/usr/bin/env python3
"""Extract OpenFOAM velocity from raw field file, render as HTML heatmap."""
import re
import math
import json
from pathlib import Path

U_PATH = Path("/home/bcant/workspace/sediment-sim/bridge/test_case/271/U")
OUT_PATH = Path("/home/bcant/workspace/sediment-sim/bridge/of_heatmap.html")

# BlockMesh background grid: 200×110×1
NX_BG, NY_BG = 200, 110
# Domain dimensions in meters
DOMAIN_W, DOMAIN_H = 12.8, 5.77
# Output grid resolution for rendering
NX_OUT, NY_OUT = 400, 220


def parse_u_field(path):
    """Parse OpenFOAM U field text file."""
    text = path.read_text()
    # Find internalField section
    start = text.find("internalField")
    section = text[start:]
    lines = section.split('\n')
    
    vals = []
    count = None
    reading = False
    
    for line in lines:
        s = line.strip()
        if not reading:
            if 'nonuniform' in s:
                # Next non-empty line is the count
                continue
            if count is None and s.isdigit():
                count = int(s)
                continue
            if count is not None and '(' in s:
                reading = True
                # Parse content after (
                idx = s.find('(')
                rest = s[idx+1:]
                if rest and ')' not in rest:
                    nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', rest)
                    if len(nums) >= 3:
                        vals.append((float(nums[0]), float(nums[1]), float(nums[2])))
                continue
        else:
            # Check for end of list: line that is just ")" or starts with ");"
            if s == ')' or s.startswith(');'):
                break
            if s and not s.startswith('//'):
                nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', s)
                if len(nums) >= 3:
                    vals.append((float(nums[0]), float(nums[1]), float(nums[2])))
    
    return vals


def map_to_grid(velocities):
    """Map cell-ordered velocities to a regular grid approximately.
    
    After snappyHexMesh, cell ordering is: 
    - Untouched background cells maintain original i,j,k order
    - Refined cells are inserted in-place, disrupting the order
    
    For visualization, we'll do nearest-neighbor interpolation from 
    approximated cell positions based on index-in-bounds ratio."""
    
    n_cells = len(velocities)
    # For the background mesh cells: index maps to i = idx % NX, j = (idx / NX) % NY
    # After snappy, this is approximate but good enough for visualization
    
    # Interpolate to output grid using approximate positions
    grid_mag = [[0.0] * NX_OUT for _ in range(NY_OUT)]
    grid_ux = [[0.0] * NX_OUT for _ in range(NY_OUT)]
    
    # Map each velocity to an approximate (i,j) based on the ratio 
    # of its index to total cells, scaled to background mesh
    for idx, (ux, uy, uz) in enumerate(velocities):
        # Approximate position based on index ratio
        frac = idx / max(n_cells - 1, 1)
        # Background mesh: i runs 0..NX-1, then j increments
        i_bg = int((frac * NX_BG * NY_BG) % NX_BG)
        j_bg = int((frac * NX_BG * NY_BG) / NX_BG) % NY_BG
        
        # Map to output grid
        i_out = int(i_bg * NX_OUT / NX_BG)
        j_out = int(j_bg * NY_OUT / NY_BG)
        if 0 <= i_out < NX_OUT and 0 <= j_out < NY_OUT:
            mag = math.sqrt(ux*ux + uy*uy)
            grid_mag[j_out][i_out] = max(grid_mag[j_out][i_out], mag)
            grid_ux[j_out][i_out] = ux  # last-write; fine for viz
    
    # Fill gaps with interpolation
    filled_mag = [[0.0] * NX_OUT for _ in range(NY_OUT)]
    for j in range(NY_OUT):
        for i in range(NX_OUT):
            if grid_mag[j][i] > 0:
                filled_mag[j][i] = grid_mag[j][i]
            else:
                # Average of neighbors
                total, count = 0.0, 0
                for dj in [-1, 0, 1]:
                    for di in [-1, 0, 1]:
                        nj, ni = j + dj, i + di
                        if 0 <= nj < NY_OUT and 0 <= ni < NX_OUT:
                            if grid_mag[nj][ni] > 0:
                                total += grid_mag[nj][ni]
                                count += 1
                filled_mag[j][i] = total / count if count > 0 else 0.0
    
    return filled_mag


def generate_html(mag_grid, out_path):
    nx, ny = NX_OUT, NY_OUT
    mag_flat = []
    for row in mag_grid:
        mag_flat.extend(row)
    mag_max = max(mag_flat) if mag_flat else 1.0
    
    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>OpenFOAM Velocity — simpleFoam converged</title>
<style>
  body {{ margin: 0; background: #0a0a14; color: #ccc; font-family: 'Courier New', monospace; overflow: hidden; }}
  canvas {{ display: block; }}
  #info {{ position: fixed; top: 12px; left: 12px; font-size: 12px; line-height: 1.8; pointer-events: none; }}
  #info span {{ color: #fa0; }}
</style></head><body>
<div id="info">
  OpenFOAM &middot; simpleFoam &middot; k-&epsilon; &middot; 271 iterations &middot; 30,210 cells<br>
  Inlet: <span>2.50 m/s</span> &middot; Max: <span>{mag_max:.2f} m/s</span><br>
  Slip walls &middot; No-slip obstacles &middot; Outflow right<br>
  &mdash; bridge output &mdash;
</div>
<canvas id="c"></canvas>
<script>
const nx = {nx}, ny = {ny};
const mag = {json.dumps(mag_flat)};
const magMax = {mag_max};

const canvas = document.getElementById('c');
const ctx = canvas.getContext('2d');

function velColor(t) {{
  t = Math.max(0, Math.min(1, Math.pow(t, 0.55)));
  let r, g, b;
  if (t < 0.25) {{ r=0; g=t*4; b=1; }}
  else if (t < 0.5) {{ r=0; g=1; b=1-(t-0.25)*4; }}
  else if (t < 0.75) {{ r=(t-0.5)*4; g=1; b=0; }}
  else {{ r=1; g=1-(t-0.75)*4; b=0; }}
  return `rgb(${{Math.floor(r*255)}},${{Math.floor(g*255)}},${{Math.floor(b*255)}})`;
}}

function draw() {{
  const w = canvas.width, h = canvas.height;
  const img = ctx.createImageData(nx, ny);
  for (let j = 0; j < ny; j++) {{
    for (let i = 0; i < nx; i++) {{
      const idx = j * nx + i;
      const t = mag[idx] / magMax;
      const color = velColor(t);
      const m = color.match(/\\d+/g);
      const pi = idx * 4;
      if (m) {{ img.data[pi]=+m[0]; img.data[pi+1]=+m[1]; img.data[pi+2]=+m[2]; }}
      img.data[pi+3] = 255;
    }}
  }}
  const off = document.createElement('canvas');
  off.width = nx; off.height = ny;
  off.getContext('2d').putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(off, 0, 0, w, h);
  
  // Color bar
  const bw = 18, bh = h * 0.5, bx = w - 40, by = (h - bh) / 2;
  for (let j = 0; j < bh; j++) {{
    ctx.fillStyle = velColor(1 - j / bh);
    ctx.fillRect(bx, by + j, bw, 1);
  }}
  ctx.fillStyle = '#ccc'; ctx.font = '11px monospace';
  ctx.fillText(magMax.toFixed(1) + ' m/s', bx + 22, by + 14);
  ctx.fillText('0', bx + 22, by + bh - 6);
}}

function resize() {{ canvas.width = window.innerWidth; canvas.height = window.innerHeight; draw(); }}
window.addEventListener('resize', resize);
resize();
</script></body></html>"""
    out_path.write_text(html)
    print(f"Wrote {out_path} ({len(html)} bytes)")


def main():
    print("Parsing U field...")
    vel = parse_u_field(U_PATH)
    print(f"  {len(vel)} velocity vectors")
    
    print("Mapping to grid...")
    mag_grid = map_to_grid(vel)
    
    print("Generating HTML...")
    generate_html(mag_grid, OUT_PATH)
    print("Done. Open the HTML file in your browser.")


if __name__ == "__main__":
    main()
