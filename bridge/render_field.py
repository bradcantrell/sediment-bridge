#!/usr/bin/env python3
"""Extract OpenFOAM velocity, interpolate to regular grid, output HTML heatmap."""
import xml.etree.ElementTree as ET
import base64
import struct
import json
import math
from pathlib import Path

VTU_PATH = Path("/home/bcant/workspace/sediment-sim/bridge/test_case/VTK/test_case_271/internal.vtu")
OUT_PATH = Path("/home/bcant/workspace/sediment-sim/bridge/of_result.html")

# Grid dimensions matching the Stam solver for comparison
NX_GRID = 320
NY_GRID = 180


def parse_vtu(vtu_path):
    tree = ET.parse(vtu_path)
    root = tree.getroot()
    piece = root.find(".//{http://www.vtk.org/VM/07082008}Piece") or root.find(".//Piece")
    n_points = int(piece.get("NumberOfPoints", 0))
    n_cells = int(piece.get("NumberOfCells", 0))

    def decode(elem, count, ncomp=3):
        raw = base64.b64decode(elem.text.strip())
        header = struct.unpack("<Q", raw[:8])[0]
        data = raw[8:8 + header]
        return struct.unpack(f"<{count * ncomp}f", data)
    
    def decode_int(elem, count, ncomp=1):
        raw = base64.b64decode(elem.text.strip())
        header = struct.unpack("<Q", raw[:8])[0]
        data = raw[8:8 + header]
        if ncomp == 1:
            return struct.unpack(f"<{count}q", data)
        else:
            return struct.unpack(f"<{count * ncomp}q", data)

    # Points (mesh vertices)
    pe = piece.find(".//{http://www.vtk.org/VM/07082008}Points") or piece.find(".//Points")
    da = pe.find(".//{http://www.vtk.org/VM/07082008}DataArray") or pe.find(".//DataArray")
    pts = decode(da, n_points)
    vertices = [(pts[i*3], pts[i*3+1], pts[i*3+2]) for i in range(n_points)]

    # Cells: get connectivity to compute centroids
    cells_elem = piece.find(".//{http://www.vtk.org/VM/07082008}Cells") or piece.find(".//Cells")
    conn_arr = None; offsets_arr = None; types_arr = None
    for arr in cells_elem.iter():
        if arr.tag.endswith("DataArray"):
            name = arr.get("Name", "")
            if name == "connectivity":
                conn_arr = arr
            elif name == "offsets":
                offsets_arr = arr
            elif name == "types":
                types_arr = arr
    
    # Parse connectivity
    total_conn = sum(1 for _ in range(n_cells))  # placeholder
    conn = decode_int(conn_arr, int(conn_arr.get("NumberOfValues", "0")), 1) if conn_arr is not None else []
    offsets = decode_int(offsets_arr, n_cells, 1) if offsets_arr is not None else []
    
    # Compute cell centers
    centers = []
    prev_offset = 0
    for off in offsets:
        # Vertices for this cell: conn[prev_offset:off]
        cell_verts = conn[prev_offset:off]
        if cell_verts:
            cx = sum(vertices[v][0] for v in cell_verts) / len(cell_verts)
            cy = sum(vertices[v][1] for v in cell_verts) / len(cell_verts)
            centers.append((cx, cy))
        else:
            centers.append((0, 0))
        prev_offset = off

    # Velocity (cell data)
    cd = piece.find(".//{http://www.vtk.org/VM/07082008}CellData") or piece.find(".//CellData")
    vel = []
    for da in cd.iter():
        if da.tag.endswith("DataArray") and da.get("Name") == "U":
            vals = decode(da, n_cells, 3)
            vel = [(vals[i*3], vals[i*3+1]) for i in range(n_cells)]
            break

    return {"centers": centers, "velocity": vel, "n_cells": n_cells, "n_points": n_points}


def interpolate_to_grid(data, nx, ny):
    """Simple nearest-neighbor interpolation to regular grid."""
    from scipy.spatial import cKDTree
    
    points = data["centers"]
    velocity = data["velocity"]
    
    # Domain bounds
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    
    # Build KD-tree for fast lookup
    tree = cKDTree(points)
    
    # Generate regular grid
    grid_u = [[0.0] * nx for _ in range(ny)]
    grid_v = [[0.0] * nx for _ in range(ny)]
    grid_mag = [[0.0] * nx for _ in range(ny)]
    
    for j in range(ny):
        y = ymax - (j + 0.5) * (ymax - ymin) / ny  # Flip y for screen coords
        for i in range(nx):
            x = xmin + (i + 0.5) * (xmax - xmin) / nx
            dist, idx = tree.query((x, y), k=1)
            u, v = velocity[idx]
            grid_u[j][i] = u
            grid_v[j][i] = v
            grid_mag[j][i] = math.sqrt(u*u + v*v)
    
    return {
        "nx": nx, "ny": ny,
        "xmin": xmin, "xmax": xmax,
        "ymin": ymin, "ymax": ymax,
        "u": grid_u, "v": grid_v, "mag": grid_mag
    }


def generate_html(grid_data, out_path):
    """Generate a self-contained HTML heatmap."""
    nx, ny = grid_data["nx"], grid_data["ny"]
    xmin, xmax = grid_data["xmin"], grid_data["xmax"]
    ymin, ymax = grid_data["ymin"], grid_data["ymax"]
    
    # Flatten velocity magnitude for JS
    mag_flat = []
    for row in grid_data["mag"]:
        mag_flat.extend(row)
    
    # Find colormap range
    mag_max = max(mag_flat)
    
    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>OpenFOAM Velocity — simpleFoam converged</title>
<style>
  body {{ margin: 0; background: #0a0a14; color: #ccc; font-family: monospace; }}
  canvas {{ display: block; }}
  #info {{ position: fixed; top: 12px; left: 12px; font-size: 12px; line-height: 1.6; }}
  #info span {{ color: #fa0; }}
</style></head><body>
<div id="info">
  OpenFOAM simpleFoam &middot; k-&epsilon; &middot; 271 iterations<br>
  Cells: 30,210 &middot; Grid: {nx}&times;{ny}<br>
  U<sub>x</sub> max: <span>{mag_max:.2f}</span> m/s &middot; Inlet: <span>2.50</span> m/s<br>
  &mdash; drag to place obstacles in live sim, bake, compare &mdash;
</div>
<canvas id="c"></canvas>
<script>
const nx = {nx}, ny = {ny};
const mag = {json.dumps(mag_flat)};
const magMax = {mag_max};
const xmin = {xmin}, xmax = {xmax}, ymin = {ymin}, ymax = {ymax};

const canvas = document.getElementById('c');
const ctx = canvas.getContext('2d');

function resize() {{
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
  draw();
}}

function velColor(t) {{
  // t in [0, 1]: blue (slow) → cyan → yellow → red (fast)
  t = Math.max(0, Math.min(1, t));
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
      // Enhance contrast
      const tc = Math.pow(t, 0.6);
      const color = velColor(tc);
      // Parse rgb string
      const m = color.match(/\\d+/g);
      const pi = (j * nx + i) * 4;
      img.data[pi] = parseInt(m[0]);
      img.data[pi+1] = parseInt(m[1]);
      img.data[pi+2] = parseInt(m[2]);
      img.data[pi+3] = 255;
    }}
  }}
  // Scale to canvas
  const offC = document.createElement('canvas');
  offC.width = nx; offC.height = ny;
  const offCtx = offC.getContext('2d');
  offCtx.putImageData(img, 0, 0);
  
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offC, 0, 0, w, h);
  
  // Overlay color bar
  const barW = 20, barH = h * 0.6;
  const barX = w - 40, barY = (h - barH) / 2;
  for (let j = 0; j < barH; j++) {{
    const t = 1 - j / barH;
    ctx.fillStyle = velColor(Math.pow(t, 0.6));
    ctx.fillRect(barX, barY + j, barW, 1);
  }}
  ctx.fillStyle = '#ccc';
  ctx.font = '11px monospace';
  ctx.fillText(magMax.toFixed(1) + ' m/s', barX + 25, barY + 12);
  ctx.fillText('0 m/s', barX + 25, barY + barH - 4);
}}

window.addEventListener('resize', resize);
resize();
</script></body></html>"""
    out_path.write_text(html)
    print(f"Wrote {out_path} ({len(html)} bytes)")


def main():
    print("Parsing VTU...")
    data = parse_vtu(VTU_PATH)
    print(f"  {data['n_cells']} cells, {len(data['centers'])} centers, {len(data['velocity'])} velocity vectors")
    
    print(f"Interpolating to {NX_GRID}×{NY_GRID} grid...")
    grid = interpolate_to_grid(data, NX_GRID, NY_GRID)
    
    print("Generating HTML...")
    generate_html(grid, OUT_PATH)
    print("Done. Open the HTML file in your browser.")


if __name__ == "__main__":
    main()
