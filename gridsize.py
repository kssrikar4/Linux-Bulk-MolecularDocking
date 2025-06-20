import numpy as np
import sys
import os

def parse_pdb(file_path):
    """Parse a PDB file and extract atom coordinates."""
    coordinates = []
    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append((x, y, z))
    return np.array(coordinates)

def calculate_box_center_and_size(coordinates):
    """Calculate the center and size of the box for AutoDock Vina."""
    min_coords = np.min(coordinates, axis=0)
    max_coords = np.max(coordinates, axis=0)
    box_size = max_coords - min_coords
    box_center = (max_coords + min_coords) / 2
    return box_center, box_size

def main():
    if len(sys.argv) != 2:
        print("Usage: python gridsize.py path/to/receptor.pdb")
        sys.exit(1)

    pdb_file = sys.argv[1]
    
    if not os.path.exists(pdb_file):
        print(f"Error: File '{pdb_file}' not found.")
        sys.exit(1)

    coordinates = parse_pdb(pdb_file)
    if coordinates.size == 0:
        print("Error: No ATOM coordinates found in the file.")
        sys.exit(1)

    box_center, box_size = calculate_box_center_and_size(coordinates)
    
    print("Grid Box Parameters for AutoDock Vina:")
    print(f"CENTER_X = {box_center[0]:.3f}")
    print(f"CENTER_Y = {box_center[1]:.3f}")
    print(f"CENTER_Z = {box_center[2]:.3f}")
    print(f"SIZE_X = {box_size[0]:.3f}")
    print(f"SIZE_Y = {box_size[1]:.3f}")
    print(f"SIZE_Z = {box_size[2]:.3f}")

if __name__ == "__main__":
    main()

