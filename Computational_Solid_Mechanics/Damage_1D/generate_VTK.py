import numpy as np

def write_vtk(Main_model, n, filename="output_timestep"):
    """
    Writes a VTK PolyData file for ParaView visualization at timestep n,
    including point data (displacements) and cell data (damage, strain, stress).

    Parameters:
      - Main_model: object with attributes:
            X: node coordinates (nNodes,)
            T: element connectivity (nEle, 2)
            U: displacement matrix (nNodes, steps)
            damage: (nEle, steps)
            strain: (nEle, steps)
            stress: (nEle, steps)
      - n: time step index
      - filename: base filename (e.g., "output_timestep_0.vtk")
    """
    # Get node coordinates and connectivity
    X = Main_model.X.flatten()         # Ensure a 1D array of node positions
    T = Main_model.T                   # Connectivity (nEle x 2)
    U_disp = Main_model.U[:, n].flatten()  # Displacements at timestep n

    n_nodes = len(X)
    n_elements = T.shape[0]

    # Get cell data for current timestep n
    damage_cell = Main_model.damage[:, n]
    strain_cell = Main_model.strain[:, n]
    stress_cell = Main_model.stress[:, n]

    with open(f"{filename}_{n}.vtk", "w") as f:
        # Header for VTK Legacy PolyData format
        f.write("# vtk DataFile Version 3.0\n")
        f.write("1D damage results in PolyData format\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        # Write points (nodes)
        f.write(f"POINTS {n_nodes} float\n")
        for x in X:
            f.write(f"{x} 0.0 0.0\n")  # Represent 1D points in 3D

        f.write("\n")
        
        # Write lines (connectivity)
        f.write(f"LINES {n_elements} {3 * n_elements}\n")
        for conn in T:
            f.write(f"2 {conn[0]} {conn[1]}\n")
        f.write("\n")
        
        # Write point data (displacements as vectors)
        f.write(f"POINT_DATA {n_nodes}\n")
        f.write("VECTORS displacement float\n")
        for u in U_disp:
            f.write(f"{u} 0.0 0.0\n")
        f.write("\n")
        
        # Write cell data for damage, strain, and stress
        f.write(f"CELL_DATA {n_elements}\n")
        
        # Damage
        f.write("SCALARS damage float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for d in damage_cell:
            f.write(f"{d}\n")
        f.write("\n")
        
        # Strain
        f.write("SCALARS strain float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for s in strain_cell:
            f.write(f"{s}\n")
        f.write("\n")
        
        # Stress
        f.write("SCALARS stress float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for s in stress_cell:
            f.write(f"{s}\n")
        f.write("\n")
        
    print(f"VTK file written to {filename}_{n}.vtk")
