#-----------------------------------------------------------------
#
#               Computational Solid Mechanics - Class 1
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#
#               FE for 1D Beam Element
#
#-----------------------------------------------------------------

import numpy as np

# Constants
# Geometry properties
A = 1          # Area of the beam (m^2)
L = 2           # Length of the beam (m)

# Elements
nEle = 10       # Number of elements
lEle = L/nEle   # Length of the elements


# Material properties
E = 210e6       # Young's Modulus
b = 2           # Non-linear parameter

# Boundary conditions
# Dirichlet BC
U0 = 0.001       # Prescribed displacement (m)

# Neumann BC 
F = 100e3       # Traction vector (N)

# Solver parameters
max_iter = 10
tol = 1e-6

class Model_part:
    def __init__(self, A, lEle, E, b, U, F, nEle):
        self.A =  A
        self.lEle = lEle
        self.E = E
        self.bterm = b
        self.pres_U = U
        self.F = F
        self.nEle = nEle            # Number of elements
        self.nNodes = nEle + 1      # Number of nodes


class Solver_params:
    def __init__(self, tol, max_iter):
        self.max_iter = max_iter        # Number of maximum iteration by the solver
        self.tol = tol              # Minimum tolerance to achieve


def compute_Kele(Main_model):
    E = Main_model.E
    A = Main_model.A
    L = Main_model.lEle

    Kele = E*A/L*np.array([[1, -1],[-1, 1]])

    return Kele


def compute_K_NLele(Main_model, uEle):
    L = Main_model.lEle
    b = Main_model.bterm

    epsilon = (uEle[1,0] - uEle[0,0])/L
    #print("strain:", epsilon)

    K_ele = E*A*(1 + 2*b*epsilon)/L*np.array([[1, -1],[-1, 1]])

    return K_ele


def assemble_Kglob(Main_model, Kele):
    nNodes = Main_model.nNodes
    nEle = Main_model.nEle

    Kglob = np.zeros((nNodes,nNodes),dtype=float)

    for i in range(0,nEle):
        Kglob[i:i+2,i:i+2] += Kele

    return Kglob


def assemble_K_NLglob(Main_model):
    nNodes = Main_model.nNodes
    nEle = Main_model.nEle

    Kglob = np.zeros((nNodes,nNodes),dtype=float)

    for i in range(0,nEle):
        node_1 = i
        node_2 = i + 2

        uEle = Main_model.U[node_1:node_2]
        #print("Local u:",uEle)

        Kele = compute_K_NLele(Main_model, uEle)
        #print(f"Kele {i}:", Kele)

        Kglob[i:i+2, i:i+2] += Kele

    return Kglob


def assemble_Fglob(Main_model):
    nNodes =  Main_model.nNodes
    F = Main_model.F

    Fglob = np.zeros((nNodes,1), dtype = float)
    Fglob[nNodes-1,0] = F


    return Fglob


def solve_NL_NewtonRaphson(Main_model, Solver, Fglob):

    max_iter = Solver.max_iter
    tol = Solver.tol

    Main_model.U = np.zeros((Main_model.nNodes,1),dtype = float)
    Main_model.U[0] = Main_model.pres_U                                               # Apply Dirichlet BC

    free_dofs = np.array([i for i in range(Main_model.nNodes) if i != 0])  # Free DOFs
    fixed_dofs = np.array([0])                                             # Fixed DOFs (Dirichlet BC)
    print("free DOFs",free_dofs,"Fixed DOFs",fixed_dofs)

    for i in range(max_iter+1):
        print(f"Non Linear iteration {i} solver")
        Kglob = assemble_K_NLglob(Main_model)
        print("Kglob", Kglob)

        R = Fglob - np.dot(Kglob, Main_model.U)
        print("Residual R", R[free_dofs])
        R[0,0] = 0.0

        if np.linalg.norm(R) < tol :  # Convergence check
            print(f"Converged in {i} iterations!")
            return Main_model.U
        
        K_uu = Kglob[free_dofs, :][:, free_dofs]
        K_uk = Kglob[free_dofs, :][:, fixed_dofs]
        
        F_u = Fglob[free_dofs] - np.dot(K_uk, Main_model.U[fixed_dofs])
        
        delta_u = np.linalg.solve(K_uu, R[free_dofs])

        Main_model.U[free_dofs] += delta_u
        Main_model.U[fixed_dofs] = Main_model.pres_U
        print(f"U {i}", Main_model.U)
        
    return Main_model.U


def solve_NewtonRaphson(Main_model, Solver, Kglob, Fglob):

    tol = Solver.tol

    U = np.zeros((Main_model.nNodes,1),dtype = float)
    U[0] = Main_model.pres_U                                               # Apply Dirichlet BC

    free_dofs = np.array([i for i in range(Main_model.nNodes) if i != 0])  # Free DOFs
    fixed_dofs = np.array([0])                                             # Fixed DOFs (Dirichlet BC)
    print("free DOFs",free_dofs,"Fixed DOFs",fixed_dofs)

    for i in range(Solver.max_iter):
        R = Fglob - np.dot(Kglob, U)
        R[0,0] = 0.0

        if np.linalg.norm(R) < tol :  # Convergence check
            print(f"Converged in {i} iterations!")
            return U
        
        K_uu = Kglob[free_dofs, :][:, free_dofs]
        K_uk = Kglob[free_dofs, :][:, fixed_dofs]
        
        F_u = Fglob[free_dofs] - np.dot(K_uk, U[fixed_dofs])
        
        delta_u = np.linalg.solve(K_uu, F_u)
        
        
        U[free_dofs] += delta_u
        
        U[fixed_dofs] = Main_model.pres_U

    return U


def main():
    Main_model = Model_part(A, lEle, E, b, U0, F, nEle)
    Solver = Solver_params(tol, max_iter)
    print("nEle:",Main_model.nEle, "nNodes:", Main_model.nNodes, "Tol:", Solver.tol, "Max Iter:", Solver.max_iter)

    #Kele = compute_Kele(Main_model)
    #print("Kele:",Kele)

    #Kglob = assemble_Kglob(Main_model, Kele)
    #print("Kglob:", Kglob)

    Fglob = assemble_Fglob(Main_model)
    print("Fglob", Fglob)

    #U = solve_NewtonRaphson(Main_model, Solver, Kglob, Fglob)
    U = solve_NL_NewtonRaphson(Main_model, Solver, Fglob)
    print("U:", U)


#-----------------------------------------------------------------


if __name__ == "__main__":
    main()