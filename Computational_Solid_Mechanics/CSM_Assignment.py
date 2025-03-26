#-----------------------------------------------------------------
#
#               Computational Solid Mechanics - Assignment damage
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#
#               FE for 1D Beam Element
#
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Model_part:
    def __init__(self):
        self.L = 1                                                  # Length of the beam (m)
        self.nEle = 5                                               # Number of elements
        self.nNodes = self.nEle + 1                                 # Number of nodes
        self.X = np.linspace(0,self.L,self.nNodes)                  # Coordinates of nodes 
        self.T = np.array([[i, i + 1] for i in range(self.nEle)])   # Connectivity of the elements

        self.fixed_dofs = np.array([0,self.nNodes-1])
        self.free_dofs = np.arange(1, self.nNodes - 1)

        # Material properties
        self.A = 0.01                           # Area of the beam (m^2)
        self.E = 70000e6                     # Young's Modulus
        self.sigma_u = 200e6                 # Ultimate stress
        self.r_0 = self.sigma_u/np.sqrt(self.E)        # Initial damage threshold
        self.beta_hs = 1.5
        self.q_inf = self.beta_hs*self.r_0
        self.A_hs = 0.1

        # Dirichlet 
        self.U_m = 0.005     # Amplitude of the prescribed displacement
        self.m = 1           # Frecuency
        self.time_span = 1   # Time span
        self.steps = 10     # Time steps
        self.vec_T = np.linspace(0,self.time_span,self.steps)
        self.U_t = self.U_m*np.sin(self.m*self.vec_T/self.time_span)

        
class Solver_params:
    def __init__(self):
        self.max_iter = 10              # Number of maximum iteration by the solver
        self.tol = 1e-8                 # Minimum tolerance to achieve

class Reference_element:
    def __init__(self, p):
        if p == 1:
            self.X_ref = np.array([-1, 1])
            self.Z_g = np.array([[1/np.sqrt(3)], [-1/np.sqrt(3)]])
            self.nZ = 2
            self.W_g = np.array([1, 1])

            self.Ni = np.hstack([(1 - self.Z_g)/2, (1 + self.Z_g)/2])
            self.Nxi = np.hstack([-1/2*np.ones((2,1)), 1/2*np.ones((2,1))])

def plot_displacement(Main_model):
    # Loop over each node

    for i in range(Main_model.nNodes):
        plt.plot(Main_model.vec_T, Main_model.U[i, :], label=f'Node {i+1}')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacement at each node over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def animate_displacement(Main_model):
    fig, ax = plt.subplots()
    ax.set_xlim(0, Main_model.L)
    ax.set_ylim(-Main_model.U_m, Main_model.U_m)
    line, = ax.plot([], [], 'bo-', label="Displacement", markersize=6)

    def update(frame):
        # Update displacement at each time step
        ax.clear()
        ax.set_xlim(0, Main_model.L)
        ax.set_ylim(-Main_model.U_m, Main_model.U_m)
        
        # Plot the nodes with displacement
        ax.plot(Main_model.X, Main_model.U[:, frame], 'bo-', label="Displacement")
        ax.set_title(f"Displacement at Time Step {frame}")
        ax.set_xlabel("Position along beam (m)")
        ax.set_ylabel("Displacement (m)")
        ax.legend()

    ani = FuncAnimation(fig, update, frames=Main_model.steps, interval=500)
    plt.show()


def Assemble_K_global(Main_model, Element):

    nNodes = Main_model.nNodes
    nEle = Main_model.nEle
    nZ = Element.nZ

    X = Main_model.X
    T = Main_model.T

    Ni = Element.Ni
    Nxi = Element.Nxi
    nZ = Element.nZ
    W_g = Element.W_g

    #Constituve model
    E = Main_model.E
    A = Main_model.A

    
    K_global = np.zeros((nNodes, nNodes))

    for ele in range(0,nEle):
        Te = T[ele,:]
        Xe = X[Te]

        h = Xe[1]- Xe[0]

        Ke = np.zeros((2,2))
        for ig in range(0, nZ):
            Nxi_ig = Nxi[ig,:]*2/h
            w_ig = W_g[ig]*h/2
            #print("Nxi_ig", Nxi_ig, "W_ig", w_ig)

            Ke += w_ig * (Nxi_ig[:, None] @ (E*A*Nxi_ig[None, :])) 
            #print("Ke", Ke)
        
        K_global[np.ix_(Te, Te)] += Ke

    return K_global

def Newton_Raphson(Main_model, Solver, Element, U_step, n):

    max_iter = Solver.max_iter
    tol = Solver.tol

    fixed_dofs = Main_model.fixed_dofs
    free_dofs = Main_model.free_dofs
    #print("Fixed DOFs:", fixed_dofs, "Free DOFs:", free_dofs)

    for i in range(max_iter + 1):
        print(f"Nonlinear iteration {i} solver")

        # Assemble global stiffness matrix
        K_global = Assemble_K_global(Main_model, Element)
        print("Global K", K_global)

        R = -np.dot(K_global, U_step)
        print("Residual",R)
        R[fixed_dofs] = 0.0

        if np.linalg.norm(R) < tol :  # Convergence check
            print(f"Converged in {i} iterations!")
            return U_step
        
        K_uu = K_global[free_dofs, :][:, free_dofs]
        
        delta_u = np.linalg.solve(K_uu, R[free_dofs])
        
        U_step[free_dofs] += delta_u
        U_step[fixed_dofs] = np.array([0.0, Main_model.U_t[n]]).reshape(-1, 1)
        print(f"U {i}", U_step)

    return U_step


def main():
    Main_model = Model_part()
    Solver = Solver_params()

    # Define element properties
    p = 1   # Order of the element
    Element = Reference_element(p)

    # Initialize displacement matrix (nNodes × time steps)
    Main_model.U = np.zeros((Main_model.nNodes, Main_model.steps), dtype=float)

    U_step = np.zeros((Main_model.nNodes, 1), dtype=float)
    for n in range(Main_model.steps):
        print(f"Time step {n}")
        U_step[-1] = Main_model.U_t[n]
        print("Imposed Ut", Main_model.U_t[n])

        # Solve nonlinear system
        U_step = Newton_Raphson(Main_model, Solver, Element, U_step, n)
        print("U_step", U_step)

        Main_model.U[:,n] = U_step.flatten()  # Store results at each time step

    plot_displacement(Main_model)
    animate_displacement(Main_model)
#-----------------------------------------------------------------


if __name__ == "__main__":
    main()