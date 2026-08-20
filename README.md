# Sediment Bridge

Real-time fluvial sediment simulation bridged to high-fidelity OpenFOAM CFD for validation. Two-phase architecture:

- **Live solver**: p5.js + Jos Stam stable fluids — real-time interactive with thousands of sediment particles
- **Physics bridge**: OpenFOAM 14 (simpleFoam via foamRun, k-epsilon) — ground-truth CFD validation

## Structure

```
├── index.html                # Live p5.js simulation (single file)
├── bridge/
│   ├── server.py             # HTTP server on :8090 (bake + status endpoints)
│   ├── case_generator.py     # OpenFOAM case generation (blockMesh, snappyHexMesh, STL obstacles)
│   ├── displace_bed.py       # Bed morphology via blockMesh point displacement
│   ├── cell_centers.py       # Velocity field extraction (foamToVTK → browser JSON)
│   ├── of_heatmap.html       # Velocity heatmap viewer
│   └── test_case/            # Example OpenFOAM case
├── README.md
├── BRIEF.md                  # Development brief (rewrite specification)
└── .gitignore
```

## Live Sim Controls

### Shapes and obstacles

| Key / Input | Action |
|-------------|--------|
| 1 | Circle obstacle |
| 2 | Rectangle obstacle |
| 3 | Freeform polygon (click vertices, click near first vertex to close) |
| Click | Place obstacle |
| Drag | Move obstacle |
| Drag handle | Resize (rect: 8 handles, circle: 4 radius handles) |
| Mouse wheel | Rotate |
| [ / ] | Grow / shrink |
| Backspace / Delete | Delete selection (obstacle, bed point, or timeline clip) |
| Escape | Deselect / cancel |

### View and state

| Key / Input | Action |
|-------------|--------|
| Space | Pause / resume |
| F | Flow field overlay |
| G | Grid overlay |
| H | Help |
| M | Bed profile editor |
| R | Reset simulation |
| S | Save PNG |
| B | Bake to OpenFOAM |

### Sidebar controls

- **Physics** — inflow velocity, settle threshold, settle probability, erode threshold, erode probability
- **Render** — trail fade, particle size, slow / fast velocity colors, deposition size, recent / oldest deposition colors, color threshold, flow scale / density / opacity, particle count
- **Toggles** — grid, flow, bed, pause, elev (spot elevations)
- **Timeline** — on/off, loop duration (30 to 600 s), mode (loop / ping-pong); choreographs obstacle appearance over a looping range
- **Bed** — 3D fishbone thumbnail with z-scale

## Bridge Pipeline

```
Live sim state → case generator → blockMesh → snappyHexMesh → simpleFoam (foamRun)
                                                                        ↓
Browser visualization ← velocity field ← cell_centers ← converged solution
```

The **bake** operation serializes the complete scene (bedform, obstacles, inflow, sediment state) and POSTs it to the bridge server at http://192.168.1.129:8090. Obstacles are exported as STL (circle, rectangle, freeform polygon) and cut with snappyHexMesh; the bed is applied via blockMesh point displacement. The returned velocity field renders alongside the real-time approximation.

## Requirements

- **Live sim** — any modern browser
- **Bridge** — OpenFOAM 14, Python 3
