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


def displace_bottom(case_dir: Path, b64_png: str, max_depth: float = 0.02):
    """Read polyMesh/points, displace bottom vertices by bed heightfield."""
    points_file = case_dir / "constant" / "polyMesh" / "points"
    if not points_file.exists():
        raise FileNotFoundError(f"No points file at {points_file}")
    
    heights = decode_png_heights(b64_png)
    in_h, in_w = len(heights), len(heights[0]) if heights else 0
    
    # Read points file
    content = points_file.read_text()
    
    # Parse the header and point list
    # Format: header ... \nN\n(\n(x y z)\n...)\n
    header_end = content.rfind('(')
    if header_end < 0:
        raise ValueError("Cannot find point list in points file")
    
    header = content[:header_end]
    body = content[header_end:]
    
    # Parse all points
    points = []
    for match in re.finditer(r'\(([-\d.e+\s]+)\)', body):
        coords = [float(x) for x in match.group(1).split()]
        if len(coords) >= 3:
            points.append(coords)
    
    if not points:
        raise ValueError("No points found")
    
    # Find domain bounds
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin = min(zs)
    
    # Displace bottom vertices (those at z ≈ zmin)
    displaced = 0
    for p in points:
        if abs(p[2] - zmin) < 1e-6:
            # Map x,y to heightfield coordinates
            hx = int((p[0] - xmin) / (xmax - xmin) * (in_w - 1)) if xmax > xmin else 0
            hy = int((ymax - p[1]) / (ymax - ymin) * (in_h - 1)) if ymax > ymin else 0  # Flip Y
            hx = max(0, min(in_w - 1, hx))
            hy = max(0, min(in_h - 1, hy))
            
            bed_h = heights[hy][hx] * max_depth
            p[2] = zmin - bed_h  # Displace downward
            displaced += 1
    
    print(f"Displaced {displaced} bottom vertices (max depth: {max_depth}m)")
    
    # Rebuild points file
    new_body = "(\n"
    for p in points:
        new_body += f"({p[0]:.10f} {p[1]:.10f} {p[2]:.10f})\n"
    new_body += ")\n"
    
    points_file.write_text(header + new_body)


if __name__ == "__main__":
    import sys
    case_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/bed_test")
    # Test with a simple dummy PNG (blank = no displacement)
    print(f"Test: {case_dir}")
    print("(Run after blockMesh to displace bottom vertices)")
