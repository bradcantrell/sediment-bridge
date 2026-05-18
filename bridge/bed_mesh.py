#!/usr/bin/env python3
"""Generate OpenFOAM polyMesh with bottom displaced by a heightfield.

Reads a depoBuffer heightfield (grayscale PNG), displaces the bottom
vertices of a structured mesh, and writes the polyMesh files directly.
"""

import struct
import base64
import math
from pathlib import Path
from io import BytesIO

# Domain dimensions (matching blockMesh)
DOMAIN_W = 12.8
DOMAIN_H = 5.77
DOMAIN_D = 0.1      # Nominal domain thickness
MAX_BED_HEIGHT = 0.02  # Max deposition depth in meters

# Mesh resolution
NX = 200
NY = 110
NZ = 1  # Single cell thick


def decode_depo_png(b64_data: str) -> list:
    """Decode a base64 PNG to a 2D grayscale heightfield.
    Returns a list of lists [NY][NX] with heights 0.0 to 1.0."""
    # Strip data URL prefix if present
    if ',' in b64_data:
        b64_data = b64_data.split(',', 1)[1]
    
    raw = base64.b64decode(b64_data)
    
    # Try reading as PNG using built-in png parser or raw RGBA
    # Since we may not have PIL, parse minimal PNG manually
    # Actually, let's just read raw RGBA pixels from the data
    # p5.js toDataURL gives a PNG — we need to decode it
    
    # Use a minimal PNG reader
    from struct import unpack
    import zlib
    
    # Check PNG signature
    if raw[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError("Not a valid PNG")
    
    # Parse IHDR chunk
    pos = 8
    width = height = 0
    pixels = None
    
    while pos < len(raw):
        length = unpack('>I', raw[pos:pos+4])[0]
        chunk_type = raw[pos+4:pos+8].decode('ascii', errors='ignore')
        
        if chunk_type == 'IHDR':
            width = unpack('>I', raw[pos+8:pos+12])[0]
            height = unpack('>I', raw[pos+12:pos+16])[0]
        
        elif chunk_type == 'IDAT':
            if pixels is None:
                pixels = raw[pos+8:pos+8+length]
            else:
                pixels += raw[pos+8:pos+8+length]
        
        elif chunk_type == 'IEND':
            break
        
        pos += 12 + length
    
    if pixels is None:
        raise ValueError("No IDAT chunks found")
    
    # Decompress
    raw_pixels = zlib.decompress(pixels)
    
    # Parse rows (each row has filter byte + RGBA pixels)
    row_size = 1 + width * 4  # filter byte + RGBA
    heights = []
    
    for y in range(height):
        row_start = y * row_size
        row_data = raw_pixels[row_start:row_start + row_size]
        row = []
        for x in range(width):
            r = row_data[1 + x*4]
            g = row_data[1 + x*4 + 1]
            b = row_data[1 + x*4 + 2]
            # Grayscale luminance
            gray = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0
            row.append(gray)
        heights.append(row)
    
    return heights


def heightfield_to_bed(heights: list, nx_out: int = NX, ny_out: int = NY) -> list:
    """Resample heightfield to mesh resolution. Returns [ny_out][nx_out] bed heights in meters."""
    in_h = len(heights)
    in_w = len(heights[0]) if in_h > 0 else 0
    
    bed = [[0.0] * nx_out for _ in range(ny_out)]
    
    for j in range(ny_out):
        # Flip Y: image y=0 is top, mesh y=0 is bottom
        src_j = int((ny_out - 1 - j) * in_h / ny_out) if in_h > 0 else 0
        src_j = max(0, min(in_h - 1, src_j))
        for i in range(nx_out):
            src_i = int(i * in_w / nx_out) if in_w > 0 else 0
            src_i = max(0, min(in_w - 1, src_i))
            bed[j][i] = heights[src_j][src_i] * MAX_BED_HEIGHT
    
    return bed


def generate_polymesh(bed: list, output_dir: Path) -> None:
    """Generate OpenFOAM polyMesh files with displaced bottom."""
    poly_dir = output_dir / "constant" / "polyMesh"
    poly_dir.mkdir(parents=True, exist_ok=True)
    
    n_points = (NX + 1) * (NY + 1) * (NZ + 1)
    n_cells = NX * NY * NZ
    n_faces = (NX + 1) * NY * NZ + NX * (NY + 1) * NZ + NX * NY * (NZ + 1)
    n_internal = (NX - 1) * NY * NZ + NX * (NY - 1) * NZ + NX * NY * (NZ - 1)
    
    dx = DOMAIN_W / NX
    dy = DOMAIN_H / NY
    dz = DOMAIN_D / NZ
    
    # ── Points ──
    with open(poly_dir / "points", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       vectorField;
    object      points;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{n_points}
(
""")
        idx = 0
        for k in range(NZ + 1):
            z = k * dz
            for j in range(NY + 1):
                y = j * dy
                for i in range(NX + 1):
                    x = i * dx
                    # Bottom vertices (k=0) get displaced by bed height
                    if k == 0 and j < NY and i < NX:
                        # Average of surrounding bed heights
                        h00 = bed[j][i] if j < NY and i < NX else 0
                        h10 = bed[j][min(i+1, NX-1)] if j < NY else 0
                        h01 = bed[min(j+1, NY-1)][i] if i < NX else 0
                        h11 = bed[min(j+1, NY-1)][min(i+1, NX-1)]
                        z_disp = -(h00 + h10 + h01 + h11) / 4.0
                    else:
                        z_disp = 0.0
                    f.write(f"({x:.6f} {y:.6f} {z + z_disp:.6f})\n")
                    idx += 1
        f.write(")\n")

    # ── Faces ──
    # Face ordering for OpenFOAM: x-normal, y-normal, z-normal
    with open(poly_dir / "faces", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       faceList;
    object      faces;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{n_faces}
(
""")
        # Helper: point index for (i,j,k)
        def pt(i, j, k):
            return k * (NX + 1) * (NY + 1) + j * (NX + 1) + i
        
        # x-normal faces (constant i)
        for k in range(NZ):
            for j in range(NY):
                for i in range(NX + 1):
                    # Face vertices: oriented so normal points FROM owner
                    if i == 0:
                        # Inlet: outward normal is -x, reverse order
                        f.write(f"4({pt(i,j,k+1)} {pt(i,j+1,k+1)} {pt(i,j+1,k)} {pt(i,j,k)})\n")
                    else:
                        # Internal or outlet: normal +x is correct
                        f.write(f"4({pt(i,j,k)} {pt(i,j+1,k)} {pt(i,j+1,k+1)} {pt(i,j,k+1)})\n")
        
        # y-normal faces (constant j)
        for k in range(NZ):
            for j in range(NY + 1):
                for i in range(NX):
                    if j == NY:
                        # Top wall: outward normal is +y, reverse order
                        f.write(f"4({pt(i,j,k+1)} {pt(i+1,j,k+1)} {pt(i+1,j,k)} {pt(i,j,k)})\n")
                    else:
                        # Internal or bottom: normal -y
                        f.write(f"4({pt(i,j,k)} {pt(i+1,j,k)} {pt(i+1,j,k+1)} {pt(i,j,k+1)})\n")
        
        # z-normal faces (constant k)
        for k in range(NZ + 1):
            for j in range(NY):
                for i in range(NX):
                    if k == 0:
                        # Bottom wall: outward normal is -z, reverse order
                        f.write(f"4({pt(i,j+1,k)} {pt(i+1,j+1,k)} {pt(i+1,j,k)} {pt(i,j,k)})\n")
                    else:
                        # Internal or top: normal +z
                        f.write(f"4({pt(i,j,k)} {pt(i+1,j,k)} {pt(i+1,j+1,k)} {pt(i,j+1,k)})\n")
        
        f.write(")\n")

    # ── Owner ──
    with open(poly_dir / "owner", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       labelList;
    object      owner;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{n_faces}
(
""")
        # x-normal faces
        # Face at i separates cell i-1 (left) and cell i (right).
        # Normal (for internal) is +x, from owner (cell i-1) to neighbour (cell i).
        # Boundary i=0: reversed, owner=cell 0. Boundary i=NX: normal +x, owner=cell NX-1.
        for k in range(NZ):
            for j in range(NY):
                for i in range(NX + 1):
                    if i == 0:
                        cell = k * NX * NY + j * NX  # cell 0
                    elif i < NX:
                        cell = k * NX * NY + j * NX + (i - 1)  # owner is left cell
                    else:
                        cell = k * NX * NY + j * NX + (NX - 1)  # outlet: last cell
                    f.write(f"{cell}\n")
        
        # y-normal faces
        # Face at j separates cell j-1 (bottom) and cell j (top).
        # Normal (for internal) is -y, from owner (cell j-1) to neighbour (cell j).
        # Boundary j=0: normal -y (outward), owner=cell 0. Boundary j=NY: reversed, owner=cell NY-1.
        for k in range(NZ):
            for j in range(NY + 1):
                for i in range(NX):
                    if j == 0:
                        cell = k * NX * NY + 0 * NX + i  # cell 0,j
                    elif j < NY:
                        cell = k * NX * NY + (j - 1) * NX + i  # owner is lower cell
                    else:
                        cell = k * NX * NY + (NY - 1) * NX + i  # top: last cell
                    f.write(f"{cell}\n")
        
        # z-normal faces: owner is cell above (higher k)
        for k in range(NZ + 1):
            for j in range(NY):
                for i in range(NX):
                    if k < NZ:
                        cell = k * NX * NY + j * NX + i
                    else:
                        cell = (k - 1) * NX * NY + j * NX + i
                    f.write(f"{cell}\n")
        
        f.write(")\n")

    # ── Neighbour ──
    with open(poly_dir / "neighbour", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       labelList;
    object      neighbour;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{n_internal}
(
""")
        # x-normal internal faces: neighbour is cell to the right (higher i)
        for k in range(NZ):
            for j in range(NY):
                for i in range(1, NX):
                    cell = k * NX * NY + j * NX + i  # cell i is the neighbour
                    f.write(f"{cell}\n")
        
        # y-normal internal faces: neighbour is cell above (higher j)
        for k in range(NZ):
            for j in range(1, NY):
                for i in range(NX):
                    cell = k * NX * NY + j * NX + i  # cell at j is the neighbour
                    f.write(f"{cell}\n")
        
        # z-normal internal faces: neighbour is cell above (higher k)
        for k in range(1, NZ):
            for j in range(NY):
                for i in range(NX):
                    cell = k * NX * NY + j * NX + i  # cell at k is the neighbour
                    f.write(f"{cell}\n")
        
        f.write(")\n")

    # ── Boundary ──
    x_faces_start = 0
    x_faces_count = (NX + 1) * NY * NZ
    y_faces_start = x_faces_start + x_faces_count
    y_faces_count = NX * (NY + 1) * NZ
    z_faces_start = y_faces_start + y_faces_count
    z_faces_count = NX * NY * (NZ + 1)
    
    with open(poly_dir / "boundary", "w") as f:
        f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    class       polyBoundaryMesh;
    object      boundary;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

5
(
    inlet
    {{
        type            patch;
        nFaces          {NY * NZ};
        startFace       {x_faces_start};
    }}
    outlet
    {{
        type            patch;
        nFaces          {NY * NZ};
        startFace       {x_faces_start + (NX) * NY * NZ};
    }}
    top
    {{
        type            wall;
        nFaces          {NX * NY};
        startFace       {z_faces_start + NZ * NX * NY};
    }}
    bottom
    {{
        type            wall;
        nFaces          {NX * NY};
        startFace       {z_faces_start};
    }}
    frontAndBack
    {{
        type            empty;
""")
        # frontAndBack: y-min and y-max faces
        f_y0_start = y_faces_start
        f_y0_count = NX * NZ
        f_yN_start = y_faces_start + NY * NX * NZ
        f_yN_count = NX * NZ
        f.write(f"        nFaces          {f_y0_count + f_yN_count};\n")
        f.write(f"        startFace       {f_y0_start};\n")
        f.write("    }\n")
        f.write(")\n")


def generate_polymesh_from_b64(b64_png: str, output_dir: Path) -> None:
    """Full pipeline: decode PNG → heightfield → polyMesh."""
    heights = decode_depo_png(b64_png)
    bed = heightfield_to_bed(heights)
    generate_polymesh(bed, output_dir)


if __name__ == "__main__":
    # Test: generate a simple sloped bed
    bed = [[0.5 * MAX_BED_HEIGHT * (1 + math.sin(i * math.pi / NX)) for i in range(NX)] for _ in range(NY)]
    test_dir = Path("/tmp/bed_test")
    generate_polymesh(bed, test_dir)
    print(f"Generated {test_dir}/constant/polyMesh/")
    print(f"Points: {(NX+1)*(NY+1)*(NZ+1)}, Cells: {NX*NY*NZ}, Faces: {(NX+1)*NY*NZ + NX*(NY+1)*NZ + NX*NY*(NZ+1)}")
