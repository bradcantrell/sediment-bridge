"""
Sediment Sim → OpenFOAM Bridge
================================
Exports live sim state, generates OpenFOAM cases, runs solver,
parses results back for comparison.

Architecture:
  Sim State (dict) → case_generator → OpenFOAM case dir
                                   → solver_runner → stdout/log
                                   → result_parser → velocity field + shear stress
"""

import json
import subprocess
import shutil
import os
import struct
import math
import tempfile
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

# ── Domain mapping ──
# p5.js canvas → OpenFOAM physical domain
# Canvas: 1280×577 px → scaled to meters for numerical stability
CANVAS_W = 1280
CANVAS_H = 577
SCALE = 0.01  # 1 px = 0.01 m → domain is 12.8m × 5.77m
DOMAIN_W = CANVAS_W * SCALE  # 12.8 m
DOMAIN_H = CANVAS_H * SCALE  # 5.77 m
DOMAIN_D = 0.1               # Thin in z (2D approximation)

# Mesh resolution
NX_BG = 200   # Background mesh cells in x
NY_BG = 110   # Background mesh cells in y
NZ_BG = 1     # Single cell in z


@dataclass
class Obstacle:
    type: str       # 'circle', 'rect', or 'poly'
    x: float        # px (centroid)
    y: float        # px (centroid)
    rotation: float = 0.0  # degrees
    w: float = 0.0  # px width (rect / circle diameter)
    h: float = 0.0  # px height (rect)
    size: float = 0.0  # legacy px (diameter or side) — fallback
    points: list = field(default_factory=list)  # polygon px points, relative to centroid


@dataclass
class SimState:
    obstacles: list = field(default_factory=list)
    inflow_velocity: float = 2.5       # m/s (scaled from grid units)
    width_px: float = CANVAS_W
    height_px: float = CANVAS_H


def px_to_m(px: float) -> float:
    """Convert pixel coordinates to meters."""
    return px * SCALE


def generate_stl_circle(cx_m: float, cy_m: float, radius_m: float,
                        z_min: float, z_max: float, n_segments: int = 32) -> str:
    """Generate ASCII STL for an extruded cylinder (circle in 2D)."""
    triangles = []
    angles = [2 * math.pi * i / n_segments for i in range(n_segments)]

    def add_triangle(v1, v2, v3, nx, ny, nz):
        triangles.append(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n"
                         f"    outer loop\n"
                         f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n"
                         f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n"
                         f"      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}\n"
                         f"    endloop\n"
                         f"  endfacet\n")

    # Side faces
    for i in range(n_segments):
        a0, a1 = angles[i], angles[(i + 1) % n_segments]
        x0, y0 = cx_m + radius_m * math.cos(a0), cy_m + radius_m * math.sin(a0)
        x1, y1 = cx_m + radius_m * math.cos(a1), cy_m + radius_m * math.sin(a1)
        nx = math.cos((a0 + a1) / 2)
        ny = math.sin((a0 + a1) / 2)
        add_triangle((x0, y0, z_min), (x1, y1, z_min), (x1, y1, z_max), nx, ny, 0)
        add_triangle((x0, y0, z_min), (x1, y1, z_max), (x0, y0, z_max), nx, ny, 0)

    # Top and bottom caps (approximate with triangle fan)
    for z, nz in [(z_min, -1.0), (z_max, 1.0)]:
        for i in range(1, n_segments - 1):
            x0, y0 = cx_m + radius_m * math.cos(angles[0]), cy_m + radius_m * math.sin(angles[0])
            xi, yi = cx_m + radius_m * math.cos(angles[i]), cy_m + radius_m * math.sin(angles[i])
            xj, yj = cx_m + radius_m * math.cos(angles[i + 1]), cy_m + radius_m * math.sin(angles[i + 1])
            if nz > 0:
                add_triangle((x0, y0, z), (xj, yj, z), (xi, yi, z), 0, 0, nz)
            else:
                add_triangle((x0, y0, z), (xi, yi, z), (xj, yj, z), 0, 0, nz)

    name = f"circle_{cx_m:.2f}_{cy_m:.2f}"
    return f"solid {name}\n" + "".join(triangles) + f"endsolid {name}\n"


def generate_stl_rect(cx_m: float, cy_m: float, half_w: float, half_h: float,
                      rotation_deg: float, z_min: float, z_max: float) -> str:
    """Generate ASCII STL for an extruded and rotated rectangle."""
    # Vertices of rectangle in local coords
    corners = [(-half_w, -half_h), (half_w, -half_h),
               (half_w, half_h), (-half_w, half_h)]
    # Rotate
    rad = math.radians(rotation_deg)
    cos_a, sin_a = math.cos(rad), math.sin(rad)
    rotated = [(cx_m + x * cos_a - y * sin_a, cy_m + x * sin_a + y * cos_a)
               for x, y in corners]

    triangles = []

    def add_triangle(v1, v2, v3, nx, ny, nz):
        triangles.append(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n"
                         f"    outer loop\n"
                         f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n"
                         f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n"
                         f"      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}\n"
                         f"    endloop\n"
                         f"  endfacet\n")

    # Side faces
    for i in range(4):
        j = (i + 1) % 4
        x0, y0 = rotated[i]
        x1, y1 = rotated[j]
        # Normal: perpendicular to edge in xy-plane
        dx, dy = x1 - x0, y1 - y0
        length = math.sqrt(dx * dx + dy * dy)
        nx, ny = -dy / length, dx / length
        add_triangle((x0, y0, z_min), (x1, y1, z_min), (x1, y1, z_max), nx, ny, 0)
        add_triangle((x0, y0, z_min), (x1, y1, z_max), (x0, y0, z_max), nx, ny, 0)

    # Top and bottom caps
    for z, nz in [(z_min, -1.0), (z_max, 1.0)]:
        if nz > 0:
            add_triangle((rotated[0][0], rotated[0][1], z),
                         (rotated[2][0], rotated[2][1], z),
                         (rotated[1][0], rotated[1][1], z), 0, 0, nz)
            add_triangle((rotated[0][0], rotated[0][1], z),
                         (rotated[3][0], rotated[3][1], z),
                         (rotated[2][0], rotated[2][1], z), 0, 0, nz)
        else:
            add_triangle((rotated[0][0], rotated[0][1], z),
                         (rotated[1][0], rotated[1][1], z),
                         (rotated[2][0], rotated[2][1], z), 0, 0, nz)
            add_triangle((rotated[0][0], rotated[0][1], z),
                         (rotated[2][0], rotated[2][1], z),
                         (rotated[3][0], rotated[3][1], z), 0, 0, nz)

    name = f"rect_{cx_m:.2f}_{cy_m:.2f}"
    return f"solid {name}\n" + "".join(triangles) + f"endsolid {name}\n"


def generate_stl_polygon(pts_m, z_min, z_max):
    """Generate ASCII STL for an extruded arbitrary polygon (pts_m = world meters, y flipped)."""
    n = len(pts_m)
    triangles = []

    def add_triangle(v1, v2, v3, nx, ny, nz):
        triangles.append(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n"
                         f"    outer loop\n"
                         f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n"
                         f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n"
                         f"      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}\n"
                         f"    endloop\n"
                         f"  endfacet\n")

    cx = sum(p[0] for p in pts_m) / n
    cy = sum(p[1] for p in pts_m) / n

    for i in range(n):
        x0, y0 = pts_m[i]
        x1, y1 = pts_m[(i + 1) % n]
        dx, dy = x1 - x0, y1 - y0
        length = math.sqrt(dx * dx + dy * dy)
        if length == 0:
            continue
        nx, ny = -dy / length, dx / length
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        if nx * (mx - cx) + ny * (my - cy) < 0:
            nx, ny = -nx, -ny
        add_triangle((x0, y0, z_min), (x1, y1, z_min), (x1, y1, z_max), nx, ny, 0)
        add_triangle((x0, y0, z_min), (x1, y1, z_max), (x0, y0, z_max), nx, ny, 0)

    for z, nz in [(z_min, -1.0), (z_max, 1.0)]:
        for i in range(n):
            x0, y0 = pts_m[i]
            x1, y1 = pts_m[(i + 1) % n]
            if nz > 0:
                add_triangle((cx, cy, z), (x0, y0, z), (x1, y1, z), 0, 0, nz)
            else:
                add_triangle((cx, cy, z), (x1, y1, z), (x0, y0, z), 0, 0, nz)

    name = f"poly_{n}"
    return f"solid {name}\n" + "".join(triangles) + f"endsolid {name}\n"


# ── OpenFOAM dictionary generators ──

def generate_blockMeshDict(state: SimState) -> str:
    """Generate blockMeshDict for a 2D-thin domain."""
    w, h, d = DOMAIN_W, DOMAIN_H, DOMAIN_D
    nx, ny, nz = NX_BG, NY_BG, NZ_BG
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
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 1;

vertices
(
    (0 0 0)              // 0
    ({w:.4f} 0 0)        // 1
    ({w:.4f} {h:.4f} 0)  // 2
    (0 {h:.4f} 0)        // 3
    (0 0 {d:.4f})        // 4
    ({w:.4f} 0 {d:.4f})  // 5
    ({w:.4f} {h:.4f} {d:.4f}) // 6
    (0 {h:.4f} {d:.4f})  // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    top
    {{
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }}
    bottom
    {{
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);

mergePatchPairs
(
);
"""


def generate_controlDict(end_time: float = 1000, write_interval: float = 100) -> str:
    """Generate controlDict for simpleFoam steady-state solver."""
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
    class       dictionary;
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time};
deltaT          1;

writeControl    runTime;
writeInterval   {write_interval};
purgeWrite      2;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
"""


def generate_fvSchemes() -> str:
    """Standard fvSchemes for steady incompressible flow."""
    return r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def generate_fvSolution() -> str:
    """fvSolution with SIMPLE settings."""
    return r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-6;
        relTol          0.01;
        smoother        GaussSeidel;
    }
    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
    k
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }
    epsilon
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
    residualControl
    {
        p               1e-4;
        U               1e-4;
        k               1e-4;
        epsilon         1e-4;
    }
}

relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; k 0.7; epsilon 0.7; }
}
"""


def generate_turbulence_properties() -> str:
    """k-epsilon turbulence model for channel flow."""
    return r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  RAS;

RAS
{
    RASModel        kEpsilon;
    turbulence      on;
    printCoeffs     on;
}
"""


def generate_momentum_transport() -> str:
    """k-epsilon turbulence model for channel flow — OpenFOAM v14 format.

    v14 consolidated turbulence dictionaries into constant/momentumTransport
    (RASModel -> model rename). Written alongside the legacy turbulenceProperties
    file for fallback compatibility.
    """
    return r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v14                                   |
|   \\  /    A nd           | Website:  www.openfoam.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      momentumTransport;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  RAS;

RAS
{
    model           kEpsilon;
    turbulence      on;
    printCoeffs     on;
}
"""


def generate_transport_properties() -> str:
    """Water transport properties."""
    return r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] 1e-06;  // water
"""


def _obstacle_boundary_entries(state: SimState, bc_type: str, bc_value: str = "") -> str:
    """Generate boundary entries for obstacle patches."""
    lines = []
    for i in range(len(state.obstacles)):
        lines.append(f"    obstacle{i}")
        lines.append("    {")
        if bc_value:
            lines.append(f"        type            {bc_type};")
            lines.append(f"        value           {bc_value};")
        else:
            lines.append(f"        type            {bc_type};")
        lines.append("    }")
    return '\n'.join(lines)


def generate_initial_U(state: SimState) -> str:
    """Generate initial/boundary velocity field with obstacle walls (no-slip)."""
    u_in = state.inflow_velocity
    obs_U = _obstacle_boundary_entries(state, "fixedValue", "uniform (0 0 0)")
    obs_block = f"\n{obs_U}\n" if obs_U else ""
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
    inlet
    {{
        type            fixedValue;
        value           uniform ({u_in} 0 0);
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    top
    {{
        type            slip;
    }}
    bottom
    {{
        type            slip;
    }}
    frontAndBack
    {{
        type            empty;
    }}{obs_block}
}}
"""


def generate_initial_p(state: SimState) -> str:
    """Pressure field — zeroGradient at obstacle walls."""
    obs_p = _obstacle_boundary_entries(state, "zeroGradient")
    obs_block = f"\n{obs_p}\n" if obs_p else ""
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
    inlet
    {{
        type            zeroGradient;
    }}
    outlet
    {{
        type            fixedValue;
        value           uniform 0;
    }}
    top
    {{
        type            slip;
    }}
    bottom
    {{
        type            slip;
    }}
    frontAndBack
    {{
        type            empty;
    }}{obs_block}
}}
"""


def generate_k(state: SimState) -> str:
    """Initial k — kqRWallFunction at obstacle walls."""
    obs_k = _obstacle_boundary_entries(state, "kqRWallFunction", "uniform 0.001")
    obs_block = f"\n{obs_k}" if obs_k else ""
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
    inlet   {{ type fixedValue; value uniform 0.001; }}
    outlet  {{ type zeroGradient; }}
    top     {{ type kqRWallFunction; value uniform 0.001; }}
    bottom  {{ type kqRWallFunction; value uniform 0.001; }}
    frontAndBack {{ type empty; }}{obs_block}
}}
"""


def generate_epsilon(state: SimState) -> str:
    """Initial epsilon — epsilonWallFunction at obstacle walls."""
    obs_eps = _obstacle_boundary_entries(state, "epsilonWallFunction", "uniform 0.0001")
    obs_block = f"\n{obs_eps}" if obs_eps else ""
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
    inlet   {{ type fixedValue; value uniform 0.0001; }}
    outlet  {{ type zeroGradient; }}
    top     {{ type epsilonWallFunction; value uniform 0.0001; }}
    bottom  {{ type epsilonWallFunction; value uniform 0.0001; }}
    frontAndBack {{ type empty; }}{obs_block}
}}
"""


def generate_nut(state: SimState) -> str:
    """Initial nut — nutkWallFunction at obstacle walls."""
    obs_nut = _obstacle_boundary_entries(state, "nutkWallFunction", "uniform 0")
    obs_block = f"\n{obs_nut}" if obs_nut else ""
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
    inlet   {{ type calculated; value uniform 0; }}
    outlet  {{ type calculated; value uniform 0; }}
    top     {{ type nutkWallFunction; value uniform 0; }}
    bottom  {{ type nutkWallFunction; value uniform 0; }}
    frontAndBack {{ type empty; }}{obs_block}
}}
"""


def generate_snappyHexMeshDict(state: SimState) -> str:
    """Generate snappyHexMeshDict for cutting obstacles."""
    stl_names = []
    for i, ob in enumerate(state.obstacles):
        stl_names.append(f'    obstacle{i}.stl {{ type triSurfaceMesh; name obstacle{i}; }}')
    stl_block = '\n'.join(stl_names) if stl_names else '    // no obstacles'

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
    class       dictionary;
    object      snappyHexMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
{stl_block}
}};

castellatedMeshControls
{{
    maxLocalCells 500000;
    maxGlobalCells 2000000;
    minRefinementCells 10;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;

    resolveFeatureAngle 30;
    planarAngle 30;
    allowFreeStandingZoneFaces true;

    features ();

    refinementRegions
    {{
    }}

    refinementSurfaces
    {{
{_refinement_surfaces(state)}
    }}

    locationInMesh (0.001 0.001 0.001);
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
    nFeatureSnapIter 10;
}}

addLayersControls
{{
    layers
    {{
    }}
    relativeSizes true;
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minFlatness 0.5;
    minVol 1e-13;
    minTetQuality 1e-9;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}}

writeFlags (scalarLevels layerSets layerFields);
mergeTolerance 1e-6;
"""


def _refinement_surfaces(state: SimState) -> str:
    lines = []
    for i in range(len(state.obstacles)):
        lines.append(f"        obstacle{i}")
        lines.append(f"        {{")
        lines.append(f"            level (2 2);")
        lines.append(f"        }}")
    return '\n'.join(lines) if lines else '        // none'


# ── Case generator ──

def create_case(state: SimState, case_dir: Path):
    """Generate a complete OpenFOAM case directory from sim state."""
    case_dir.mkdir(parents=True, exist_ok=True)

    # System directory
    sys_dir = case_dir / "system"
    sys_dir.mkdir(exist_ok=True)
    (sys_dir / "controlDict").write_text(generate_controlDict())
    (sys_dir / "fvSchemes").write_text(generate_fvSchemes())
    (sys_dir / "fvSolution").write_text(generate_fvSolution())
    (sys_dir / "blockMeshDict").write_text(generate_blockMeshDict(state))
    (sys_dir / "snappyHexMeshDict").write_text(generate_snappyHexMeshDict(state))
    # surfaceFeatureExtractDict — required by snappy even if empty
    (sys_dir / "surfaceFeatureExtractDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object surfaceFeatureExtractDict; }\n\n"
        "// No feature edges extracted\n"
    )

    # Constant directory
    const_dir = case_dir / "constant"
    const_dir.mkdir(exist_ok=True)
    (const_dir / "transportProperties").write_text(generate_transport_properties())
    (const_dir / "turbulenceProperties").write_text(generate_turbulence_properties())
    # OpenFOAM v14 format (constant/momentumTransport) — written alongside the
    # legacy file for fallback compatibility (v14 reads momentumTransport first).
    (const_dir / "momentumTransport").write_text(generate_momentum_transport())

    # TriSurface directory for STL files
    tri_dir = const_dir / "triSurface"
    tri_dir.mkdir(exist_ok=True)
    z_min, z_max = -DOMAIN_D / 2, DOMAIN_D / 2
    for i, ob in enumerate(state.obstacles):
        cx_m = px_to_m(ob.x)
        cy_m = px_to_m(DOMAIN_H / SCALE - ob.y)  # Flip y: canvas y=0 → top, OpenFOAM y=0 → bottom
        if ob.type == "circle":
            r_m = px_to_m((ob.w or ob.size) / 2)
            stl_content = generate_stl_circle(cx_m, cy_m, r_m, z_min, z_max)
        elif ob.type == "poly":
            a = math.radians(ob.rotation or 0.0)
            cosA, sinA = math.cos(a), math.sin(a)
            pts_m = []
            for px, py in ob.points:
                wx = ob.x + px * cosA - py * sinA
                wy = ob.y + px * sinA + py * cosA
                pts_m.append((px_to_m(wx), px_to_m(DOMAIN_H / SCALE - wy)))
            stl_content = generate_stl_polygon(pts_m, z_min, z_max)
        else:
            hw = px_to_m((ob.w or ob.size) / 2)
            hh = px_to_m((ob.h or ob.size) / 2)
            stl_content = generate_stl_rect(cx_m, cy_m, hw, hh, ob.rotation, z_min, z_max)
        (tri_dir / f"obstacle{i}.stl").write_text(stl_content)

    # Time directory 0
    t0_dir = case_dir / "0"
    t0_dir.mkdir(exist_ok=True)
    (t0_dir / "U").write_text(generate_initial_U(state))
    (t0_dir / "p").write_text(generate_initial_p(state))
    (t0_dir / "k").write_text(generate_k(state))
    (t0_dir / "epsilon").write_text(generate_epsilon(state))
    (t0_dir / "nut").write_text(generate_nut(state))


# ── Solver runner ──

def run_case(case_dir: Path) -> dict:
    """Run blockMesh → snappyHexMesh → simpleFoam in sequence.
    Returns dict with success status and output."""
    results = {"success": False, "steps": {}}

    # Need to source OpenFOAM environment first
    of_prefix = "export ParaView_TYPE=none && source /opt/openfoam14/etc/bashrc && "

    steps = [
        ("blockMesh", f"{of_prefix} blockMesh -case {case_dir}"),
        ("snappyHexMesh", f"{of_prefix} snappyHexMesh -overwrite -case {case_dir}"),
        ("simpleFoam", f"{of_prefix} simpleFoam -case {case_dir}"),
    ]

    for name, cmd in steps:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
        results["steps"][name] = {
            "exit_code": proc.returncode,
            "stderr": proc.stderr[-500:],  # Last 500 chars
        }
        if proc.returncode != 0:
            results["error"] = f"{name} failed with exit code {proc.returncode}"
            return results

    results["success"] = True
    return results


# ── Result parser ──

def parse_velocity_field(case_dir: Path, time: str = "1000") -> Optional[dict]:
    """Parse the U field from an OpenFOAM time directory.
    Returns dict with velocity data for comparison with Stam solver."""
    u_file = case_dir / time / "U"
    if not u_file.exists():
        # Try the latest time directory
        time_dirs = sorted(
            [d for d in case_dir.iterdir() if d.is_dir() and d.name.replace('.', '').isdigit()],
            key=lambda d: float(d.name)
        )
        if not time_dirs:
            return None
        u_file = time_dirs[-1] / "U"
        if not u_file.exists():
            return None

    content = u_file.read_text()
    # Parse the OpenFOAM field format
    # We want: dimensions, internalField values, and boundary values
    result = {"dimensions": None, "internal_field": [], "boundaries": {}}

    lines = content.split('\n')
    i = 0
    # Skip header
    while i < len(lines) and not lines[i].startswith("dimensions"):
        i += 1
    if i < len(lines):
        result["dimensions"] = lines[i].strip()
    i += 1

    # Find internalField
    while i < len(lines) and "internalField" not in lines[i]:
        i += 1
    if i < len(lines):
        field_line = lines[i]
        if "nonuniform" in field_line:
            # Nonuniform list: next line is count, then values
            i += 1
            while i < len(lines) and not lines[i].strip().isdigit():
                i += 1
            count = int(lines[i].strip()) if i < len(lines) else 0
            i += 1
            # Skip '('
            while i < len(lines) and '(' not in lines[i]:
                i += 1
            i += 1
            # Read values until ')'
            parsed = 0
            while i < len(lines) and parsed < count:
                line = lines[i].strip()
                if ')' in line:
                    break
                # Values are like "(2.5 0 0)" or just "2.5 0 0"
                if line and not line.startswith('//'):
                    parts = line.replace('(', '').replace(')', '').split()
                    if len(parts) >= 3:
                        try:
                            result["internal_field"].append({
                                "u": float(parts[0]),
                                "v": float(parts[1]),
                                "w": float(parts[2])
                            })
                            parsed += 1
                        except ValueError:
                            pass
                i += 1
        elif "uniform" in field_line:
            # Uniform value
            parts = field_line.split()
            val_idx = parts.index("uniform") + 1
            val_str = ' '.join(parts[val_idx:]).replace('(', '').replace(')', '')
            vals = [float(x) for x in val_str.split()]
            result["internal_field"] = {"uniform": vals}

    return result


def extract_midplane_velocity(case_dir: Path) -> Optional[list]:
    """Extract 2D velocity field at midplane for comparison with Stam solver.
    Returns a list of dicts with x, y, u, v."""
    # For a 1-cell-thick domain, all cells are at the midplane.
    # We just need the cell center positions.
    # Use foamToVTK to convert, or parse cell centres from the solver output.
    # For now, parse internal field directly.
    field_data = parse_velocity_field(case_dir)
    if not field_data or not field_data["internal_field"]:
        return None

    internal = field_data["internal_field"]
    if isinstance(internal, dict) and "uniform" in internal:
        # Uniform field — not useful
        return None

    # Need cell centres — run postProcess to get them
    # For a first pass, we can approximate: blockMesh with simpleGrading
    # gives uniform cells. We know NX_BG × NY_BG × NZ_BG cells.
    # After snappyHexMesh, cell count changes. We'd need to parse.
    return internal  # Returns raw velocity values for now


# ── Main export interface ──

def bridge_run(state: SimState, case_dir: Optional[Path] = None) -> dict:
    """Main entry point: export state → run OpenFOAM → return results.

    Args:
        state: SimState with obstacles and inflow parameters
        case_dir: Optional path for the case (uses temp dir if None)

    Returns:
        dict with success, velocity_field, error info
    """
    if case_dir is None:
        case_dir = Path(tempfile.mkdtemp(prefix="sediment_of_"))

    # Generate case
    create_case(state, case_dir)

    # Run solver
    results = run_case(case_dir)

    if results["success"]:
        # Parse velocity field
        vel = parse_velocity_field(case_dir)
        results["velocity_field"] = vel

    results["case_dir"] = str(case_dir)
    return results


# ── CLI ──

if __name__ == "__main__":
    # Test with hardcoded state matching the current sim's 3 obstacles
    ghostSize = 40
    state = SimState(
        obstacles=[
            Obstacle(type="circle", x=256, y=202, size=ghostSize),
            Obstacle(type="square", x=448, y=346, size=ghostSize, rotation=30),
            Obstacle(type="square", x=704, y=231, size=ghostSize, rotation=-15),
        ],
        inflow_velocity=2.5
    )
    test_dir = Path("/home/bcant/workspace/sediment-sim/bridge/test_case")
    print(f"Generating test case in {test_dir}...")
    create_case(state, test_dir)
    print("Done. Case directory created.")
    print(f"\nTo run: cd {test_dir} && blockMesh && snappyHexMesh -overwrite && simpleFoam")
