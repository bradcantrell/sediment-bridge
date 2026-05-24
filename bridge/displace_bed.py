#!/usr/bin/env python3
"""Displace bottom vertices of an existing OpenFOAM polyMesh using a bed heightfield.

Much simpler than generating the full mesh — blockMesh handles the valid topology,
we just move the z-coordinate of bottom-face vertices.
"""

import struct
import base64
import zlib
import re
from pathlib import Path


def decode_png_heights(b64_data: str) -> list:
    """Decode a base64 PNG to 2D grayscale heights [0..1]."""
    if ',' in b64_data:
        b64_data = b64_data.split(',', 1)[1]
    raw = base64.b64decode(b64_data)
    
    if raw[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError("Not a valid PNG")
    
    pos = 8
    width = height = 0
    pixels = None
    
    while pos < len(raw):
        length = struct.unpack('>I', raw[pos:pos+4])[0]
        chunk_type = raw[pos+4:pos+8].decode('ascii', errors='ignore')
        if chunk_type == 'IHDR':
            width = struct.unpack('>I', raw[pos+8:pos+12])[0]
            height = struct.unpack('>I', raw[pos+12:pos+16])[0]
        elif chunk_type == 'IDAT':
            pixels = (pixels or b'') + raw[pos+8:pos+8+length]
        elif chunk_type == 'IEND':
            break
        pos += 12 + length
    
    if pixels is None:
        raise ValueError("No IDAT chunks")
    
    raw_pixels = zlib.decompress(pixels)
    row_size = 1 + width * 4
    heights = []
    
    for y in range(height):
        row_data = raw_pixels[y * row_size:(y + 1) * row_size]
        row = []
        for x in range(width):
            r = row_data[1 + x*4]
            g = row_data[1 + x*4 + 1]
            b = row_data[1 + x*4 + 2]
            row.append((0.299 * r + 0.587 * g + 0.114 * b) / 255.0)
        heights.append(row)
    
    return heights


def profile_to_heights(bed_profile: list, nx: int, ny: int, domain_w: float, domain_h: float,
                       max_depth: float = 0.02, heights_png: list = None) -> list:
    """Convert a bed profile (list of {x, z} px points) to a 2D heightfield [ny][nx].
    The profile defines z(x) in pixel space; we interpolate to the full grid.
    PNG deposition data is added on top if provided.
    Returns list of lists [ny][nx] with heights in meters.
    """
    points = bed_profile  # [{x: px, z: px}]
    if len(points) < 2:
        # No profile — use PNG only or flat
        if heights_png:
            return [[heights_png[y][x] * max_depth for x in range(nx)] for y in range(ny)]
        return [[0.0] * nx for _ in range(ny)]
    
    # Sort by x
    sorted_pts = sorted(points, key=lambda p: p['x'])
    
    # Build interpolated profile z(x) at mesh resolution
    profile_z = [0.0] * nx
    for i in range(nx):
        x_px = i * domain_w / nx  # mesh x in meters → back to px for profile lookup
        # Find bracketing points
        if x_px <= sorted_pts[0]['x']:
            z_px = sorted_pts[0]['z']
        elif x_px >= sorted_pts[-1]['x']:
            z_px = sorted_pts[-1]['z']
        else:
            # Linear interpolation
            for k in range(len(sorted_pts) - 1):
                if sorted_pts[k]['x'] <= x_px <= sorted_pts[k + 1]['x']:
                    t = (x_px - sorted_pts[k]['x']) / (sorted_pts[k + 1]['x'] - sorted_pts[k]['x'])
                    z_px = sorted_pts[k]['z'] + t * (sorted_pts[k + 1]['z'] - sorted_pts[k]['z'])
                    break
        # Convert px to meters
        profile_z[i] = z_px * (domain_h / 577.0)  # Scale from px to meters using domain height
    
    # Build full heightfield
    bed = [[0.0] * nx for _ in range(ny)]
    for j in range(ny):
        for i in range(nx):
            h = profile_z[i]  # Base profile (varies only with x)
            if heights_png:
                h += heights_png[j][i] * max_depth  # Add deposition detail
            bed[j][i] = h
    
    return bed


def displace_bottom(case_dir: Path, b64_png: str, max_depth: float = 0.02,
                    bed_profile: list = None):
    """Read polyMesh/points, displace bottom vertices by bed heightfield.
    
    Args:
        case_dir: OpenFOAM case directory (with constant/polyMesh/points)
        b64_png: Deposition buffer as base64 PNG (can be empty string)
        max_depth: Maximum deposition depth in meters
        bed_profile: List of {x, z} pixel coords for user-drawn bed profile
    """
    points_file = case_dir / "constant" / "polyMesh" / "points"
    if not points_file.exists():
        raise FileNotFoundError(f"No points file at {points_file}")
    
    # Decode PNG heights if provided
    png_heights = None
    if b64_png and len(b64_png) > 100:
        try:
            png_heights = decode_png_heights(b64_png)
        except Exception:
            pass  # Use flat if PNG decode fails
    
    in_h = len(png_heights) if png_heights else 0
    in_w = len(png_heights[0]) if png_heights else 0
    
    # Read points file
    content = points_file.read_text()
    
    # Find the opening parenthesis of the point list
    # Format: ...header... \nN\n(\n(x y z)\n...)\n
    # Find the line that is just "(" after the count
    lines = content.split('\n')
    open_idx = -1
    for i, line in enumerate(lines):
        s = line.strip()
        if s == '(' and i > 0 and lines[i-1].strip().isdigit():
            open_idx = i
            break
    
    if open_idx < 0:
        # Fallback: find first '(' after "points" keyword
        for i, line in enumerate(lines):
            if line.strip() == '(':
                open_idx = i
                break
    
    if open_idx < 0:
        raise ValueError("Cannot find opening parenthesis in points file")
    
    # Header = everything up to and including the opening '(' line
    header_lines = lines[:open_idx + 1]
    body_lines = lines[open_idx + 1:]
    
    # Parse all points from body
    points = []
    for line in body_lines:
        s = line.strip()
        if s == ')':
            break
        nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', s)
        if len(nums) >= 3:
            points.append([float(nums[0]), float(nums[1]), float(nums[2])])
    
    if not points:
        raise ValueError("No points found")
    
    # Find domain bounds
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin = min(zs)
    
    # Number of grid cells in x and y (from the structured mesh)
    # Approximate from the point count: sqrt of total points ≈ (nx+1)*(ny+1)
    # For our standard 200×110 mesh: nx=200, ny=110
    nx_mesh = 200
    ny_mesh = 110
    
    # Build heightfield from profile + PNG
    if bed_profile and len(bed_profile) >= 2:
        heights_grid = profile_to_heights(bed_profile, nx_mesh, ny_mesh,
                                          xmax - xmin, ymax - ymin,
                                          max_depth, png_heights)
    elif png_heights:
        # PNG only, resample to mesh resolution
        heights_grid = [[0.0] * nx_mesh for _ in range(ny_mesh)]
        for j in range(ny_mesh):
            for i in range(nx_mesh):
                if in_w > 0 and in_h > 0:
                    pi = int(i * in_w / nx_mesh)
                    pj = int(j * in_h / ny_mesh)
                    pi = max(0, min(in_w - 1, pi))
                    pj = max(0, min(in_h - 1, pj))
                    heights_grid[j][i] = png_heights[pj][pi] * max_depth
    else:
        heights_grid = [[0.0] * nx_mesh for _ in range(ny_mesh)]
    
    # Displace bottom vertices (those at z ≈ zmin)
    displaced = 0
    for p in points:
        if abs(p[2] - zmin) < 1e-6:
            # Map x,y to mesh indices
            mi = int((p[0] - xmin) / (xmax - xmin) * (nx_mesh - 1)) if xmax > xmin else 0
            mj = int((p[1] - ymin) / (ymax - ymin) * (ny_mesh - 1)) if ymax > ymin else 0
            mi = max(0, min(nx_mesh - 1, mi))
            mj = max(0, min(ny_mesh - 1, mj))
            
            bed_h = heights_grid[mj][mi]
            p[2] = zmin - bed_h  # Displace downward
            displaced += 1
    
    print(f"Displaced {displaced} bottom vertices (max depth: {max_depth}m)")
    
    # Rebuild points file
    result_lines = list(header_lines)
    # header_lines already includes the opening '(' — don't add another
    for p in points:
        result_lines.append(f"({p[0]:.10f} {p[1]:.10f} {p[2]:.10f})")
    result_lines.append(")")
    result_lines.append("")
    
    points_file.write_text('\n'.join(result_lines))


if __name__ == "__main__":
    import sys
    case_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/bed_test")
    # Test with a simple dummy PNG (blank = no displacement)
    print(f"Test: {case_dir}")
    print("(Run after blockMesh to displace bottom vertices)")
