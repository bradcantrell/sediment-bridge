# Sediment Simulation Rewrite — Development Brief

## Current State Analysis

The existing sediment-sim at `~/workspace/sediment-sim/index.html` (1148 lines, p5.js + Stam fluid solver) has several structural issues in the particle simulation:

### Critical Issues
1. **No particle inertia**: `updateParticles()` line 444 directly assigns `p.vx = vel.vx; p.vy = vel.vy;` — particles instantly adopt flow velocity. No mass, no momentum, no drag. Particles should gradually accelerate toward flow velocity with a response coefficient.
2. **Position noise instead of velocity noise**: Lines 461-462 add Perlin noise directly to position (`p.x += noise(...)`), causing unnatural jittering. Noise should perturb velocity (turbulent diffusion), not position.
3. **DT multiplier inconsistency**: Line 459 `p.x += p.vx * DT * 4` — the DT*4 factor is arbitrary. Particle advection should match the Stam solver's DT convention.
4. **Floor-based cell indexing**: `getVelocityAt` uses `floor(px / cellW)` which can produce off-by-one indexing errors at grid boundaries when `cellW = width/NX` rather than `width/(NX+2)`.

### Secondary Issues
5. **No particle age culling**: `deposited` array grows unboundedly. Particles older than some threshold that never deposit should be culled.
6. **Simple obstacle collision**: `isInsideObstacle` just pushes particles randomly instead of proper reflection/sliding.
7. **Constant inflow velocity**: Inflow is uniform across the left boundary with small noise. Should support variable inflow profiles.
8. **No restart capability**: Cannot reset/reinitialize particles without page reload.

## Rewrite Specification

### Architecture
Rewrite as a **clean, single-file p5.js application** at `~/workspace/sediment-sim/index.html`. Preserve the existing UI (sidebar controls, HUD, keyboard shortcuts, bed profile editor, bake button, fishbone thumbnail). The rewrite should be a complete replacement of the current file.

### Technical Requirements
1. **Stam solver**: Keep the existing 160×90 Stam grid. Fix `cellW = width / (NX + 2)`. Ensure boundary conditions are correctly applied.
2. **Particle inertia model**: Particles have mass=1. Apply flow force: `dv = (flowVel - particleVel) * responseCoeff * DT`. Response coefficient ~0.3-0.5. Add optional drag term.
3. **Turbulent diffusion**: Add velocity noise proportional to flow speed: `vx += (noise(...) - 0.5) * turbulenceStrength * speed`. No position noise.
4. **Particle advection**: Use consistent DT: `p.x += p.vx * DT; p.y += p.vy * DT;` (no arbitrary multiplier). Use RK2 or at minimum velocity-at-midpoint for accuracy.
5. **Particle culling**: Remove particles that travel beyond canvas bounds OR live longer than MAX_AGE (30 seconds) without depositing.
6. **Obstacle collision**: Push particles to nearest non-obstacle position, reflect velocity component normal to surface.
7. **Restart**: Press 'R' resets all particles and simulation state.
8. **Performance**: Maintain 30+ FPS with 3000 particles. Use typed arrays where beneficial. Keep GAUSS_SEIDEL_ITERS at 30.

### What to Preserve
- Sidebar with all existing controls (physics sliders, render sliders, toggles, color pickers)
- HUD overlay (sim-time, particles, deposited, obstacles)
- Keyboard shortcuts (F=flow, G=grid, M=bed, Space=pause, R=restart, S=save, B=bake)
- Bed profile editor with longitudinal + transverse points
- Fishbone 3D thumbnail with rotation and z-scale
- Trail buffer (persistence-of-vision)
- Flow vector overlay
- Bake-to-OpenFOAM bridge (POST to :8090)
- Dark theme aesthetic

### Code Quality
- Clean separation of concerns: solver, particles, rendering, UI
- No global mutable state beyond what's necessary
- Clear variable naming
- Comments explaining the physics, not the code
