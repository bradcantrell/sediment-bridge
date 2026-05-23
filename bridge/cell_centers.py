def compute_cell_centers(poly_mesh_dir):
    """Compute cell centers from OpenFOAM polyMesh files.
    Returns list of (x, y) tuples in meters.
    """
    import re
    
    # Read points
    points_text = (poly_mesh_dir / "points").read_text()
    points = []
    in_list = False
    for line in points_text.split('\n'):
        s = line.strip()
        if s == '(':
            in_list = True
            continue
        if in_list:
            if s == ')':
                break
            nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', s)
            if len(nums) >= 3:
                points.append((float(nums[0]), float(nums[1]), float(nums[2])))
    
    # Read faces
    faces_text = (poly_mesh_dir / "faces").read_text()
    faces = []
    in_list = False
    current_face = []
    for line in faces_text.split('\n'):
        s = line.strip()
        if s == '(':
            in_list = True
            continue
        if in_list:
            if s == ')':
                break
            # Format: "4(0 1 2 3)"
            m = re.match(r'\d+\((.*)\)', s)
            if m:
                verts = [int(x) for x in m.group(1).split()]
                faces.append(verts)
    
    # Read owner
    owner_text = (poly_mesh_dir / "owner").read_text()
    owners = []
    for line in owner_text.split('\n'):
        s = line.strip()
        if s.isdigit() or (s.startswith('-') and s[1:].isdigit()):
            owners.append(int(s))
    
    # Read neighbour (only internal faces have neighbours)
    neighbour_file = poly_mesh_dir / "neighbour"
    neighbours = []
    if neighbour_file.exists():
        neighbour_text = neighbour_file.read_text()
        for line in neighbour_text.split('\n'):
            s = line.strip()
            if s.isdigit() or (s.startswith('-') and s[1:].isdigit()):
                neighbours.append(int(s))
    
    # Determine number of cells: count unique owners (not max+1, cells may have gaps)
    n_cells = len(set(owners))
    # Also include cells referenced in neighbours that might not be in owners
    if neighbours:
        all_cells = set(owners) | set(neighbours)
        n_cells = len(all_cells)
    else:
        n_cells = len(set(owners))
    
    # For each cell, collect face centers (use dict since cell indices may have gaps)
    from collections import defaultdict
    cell_face_centers = defaultdict(list)
    
    for face_idx, vert_indices in enumerate(faces):
        if not vert_indices:
            continue
        # Face center
        cx = sum(points[v][0] for v in vert_indices) / len(vert_indices)
        cy = sum(points[v][1] for v in vert_indices) / len(vert_indices)
        center = (cx, cy)
        
        # Owner cell
        if face_idx < len(owners):
            owner = owners[face_idx]
            cell_face_centers[owner].append(center)
        
        # Neighbour cell (internal faces only)
        if face_idx < len(neighbours):
            neighbour = neighbours[face_idx]
            cell_face_centers[neighbour].append(center)
    
    # Cell center = average of face centers, ordered by cell index
    all_cells = sorted(cell_face_centers.keys())
    cell_centers = []
    for cell_idx in all_cells:
        face_centers = cell_face_centers[cell_idx]
        cx = sum(fc[0] for fc in face_centers) / len(face_centers)
        cy = sum(fc[1] for fc in face_centers) / len(face_centers)
        cell_centers.append((cx, cy))
    
    return cell_centers
