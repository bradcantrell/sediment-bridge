#!/usr/bin/env python3
"""Extract OpenFOAM velocity field from VTK output for comparison with Stam solver."""
import xml.etree.ElementTree as ET
import base64
import struct
import json
from pathlib import Path

VTU_PATH = Path("/home/bcant/workspace/sediment-sim/bridge/test_case/VTK/test_case_271/internal.vtu")


def parse_vtu(vtu_path: Path) -> dict:
    """Parse a VTU unstructured grid file and extract points + velocity field."""
    tree = ET.parse(vtu_path)
    root = tree.getroot()

    # Find the Piece element
    piece = root.find(".//{http://www.vtk.org/VM/07082008}Piece")
    if piece is None:
        piece = root.find(".//Piece")
    
    n_points = int(piece.get("NumberOfPoints", 0))
    n_cells = int(piece.get("NumberOfCells", 0))

    # Find PointData and CellData
    ns = "{http://www.vtk.org/VM/07082008}"

    def decode_data_array(elem, expected_count, n_components=3):
        """Decode a VTK binary DataArray with UInt64 header."""
        raw = base64.b64decode(elem.text.strip())
        # First 8 bytes: uint64 header with byte count
        header_size = struct.unpack("<Q", raw[:8])[0]
        data = raw[8:8 + header_size]
        fmt = f"<{expected_count * n_components}f"
        vals = struct.unpack(fmt, data)
        return vals
    
    # Points (cell centers for the internal mesh)
    points_elem = piece.find(f".//{ns}Points")
    if points_elem is None:
        points_elem = piece.find(".//Points")
    
    points = []
    if points_elem is not None:
        data_arr = points_elem.find(f".//{ns}DataArray[@Name='Points']")
        if data_arr is None:
            data_arr = points_elem.find(".//DataArray[@Name='Points']")
        if data_arr is not None:
            pts = decode_data_array(data_arr, n_points, 3)
            points = [(pts[i*3], pts[i*3+1], pts[i*3+2]) for i in range(n_points)]

    # Cell data (velocity is cell-centered for OpenFOAM internal fields)
    cell_data = piece.find(f".//{ns}CellData")
    if cell_data is None:
        cell_data = piece.find(".//CellData")
    
    velocity = []
    if cell_data is not None:
        for da in cell_data.iter():
            if da.tag.endswith("DataArray") and da.get("Name") == "U":
                n_comp = int(da.get("NumberOfComponents", 3))
                vals = decode_data_array(da, n_cells, n_comp)
                velocity = [(vals[i*3], vals[i*3+1], vals[i*3+2]) for i in range(n_cells)]
                break

    return {
        "n_cells": n_cells,
        "n_points": n_points,
        "velocity": velocity,
        "points": points,
    }


def velocity_stats(velocity: list) -> dict:
    """Compute statistics from velocity field."""
    if not velocity:
        return {}
    ux = [v[0] for v in velocity]
    uy = [v[1] for v in velocity]
    return {
        "n_cells": len(velocity),
        "ux_min": min(ux),
        "ux_max": max(ux),
        "ux_mean": sum(ux) / len(ux),
        "uy_max_abs": max(abs(v) for v in uy),
        "reverse_flow": sum(1 for v in ux if v < 0),
        "stagnant": sum(1 for v in ux if abs(v) < 0.01),
        "below_0_5": sum(1 for v in ux if v < 0.5),
    }


def main():
    data = parse_vtu(VTU_PATH)
    stats = velocity_stats(data["velocity"])
    
    # Sample: first 5, middle 5, last 5 velocities
    v = data["velocity"]
    n = len(v)
    samples = {
        "first_5": v[:5],
        "mid_5": v[n//2-2:n//2+3],
        "last_5": v[-5:],
    }
    
    output = {"stats": stats, "samples": samples}
    print(json.dumps(output, indent=2, default=str))


if __name__ == "__main__":
    main()
