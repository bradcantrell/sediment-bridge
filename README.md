# Sediment Sim

Interactive sediment transport simulation with two-phase architecture:

- **Live solver**: p5.js + Jos Stam stable fluids — real-time interactive at 60 FPS
- **Physics bridge**: OpenFOAM (simpleFoam, k-epsilon) — ground-truth CFD validation

## Structure

```
├── artworks/sediment-sim/
│   └── index.html          # Live p5.js simulation
├── bridge/
│   ├── case_generator.py   # OpenFOAM case generation from sim state
│   ├── render_heatmap.py   # Velocity field visualization
│   └── test_case/          # Example OpenFOAM case
└── .gitignore
```

## Live Sim Controls

| Key | Action |
|-----|--------|
| 1/2 | Circle / Square obstacle |
| Click | Place obstacle |
| Drag | Move obstacle |
| Backspace | Delete selected |
| Mouse wheel | Rotate |
| Space | Pause |
| F | Flow field overlay |
| G | Grid overlay |
| R | Reset |
| S | Save PNG |

## Bridge Pipeline

```
Live sim state → Python case generator → blockMesh → snappyHexMesh → simpleFoam
                                                                       ↓
Browser visualization ← velocity field ← foamToVTK ← converged solution
```

## Requirements

- **Live sim**: Any modern browser
- **Bridge**: OpenFOAM v1912+, Python 3, scipy
