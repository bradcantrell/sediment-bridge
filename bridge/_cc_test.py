import sys
sys.path.insert(0, '/home/bcant/workspace/sediment-sim/bridge')
from cell_centers import compute_cell_centers
from pathlib import Path

centers = compute_cell_centers(Path('/home/bcant/workspace/sediment-sim/bridge/bake_case/constant/polyMesh'))
print(f'Cell centers: {len(centers)}')
if centers:
    print(f'First: {centers[0]}')
    xs = [c[0] for c in centers]
    ys = [c[1] for c in centers]
    print(f'X: {min(xs):.3f} - {max(xs):.3f}')
    print(f'Y: {min(ys):.3f} - {max(ys):.3f}')
print(f'\nVelocity vectors expected: 27716')
print(f'Match: {len(centers) == 27716}')
