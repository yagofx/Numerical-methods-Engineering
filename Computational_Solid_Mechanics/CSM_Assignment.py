#-----------------------------------------------------------------
#
#               Computational Solid Mechanics - Assignment damage
#                   
#               Authors: Gabriel Ayu, David Conteras, Yago Trias
#               Co-authors: Matteo Giaccomini (el pro)
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#
#               FE for 1D Beam damage model
#
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Model_part:
    def __init__(self):
        self.L = 1                                                  # Length of the beam (m)
        self.nEle = 6                                               # Number of elements
        self.nNodes = self.nEle + 1                                 # Number of nodes
        self.X = np.linspace(0,self.L,self.nNodes)                  # Coordinates of nodes 
        self.T = np.array([[i, i + 1] for i in range(self.nEle)])   # Connectivity of the elements

        self.fixed_dofs = np.array([0,self.nNodes-1])
        self.free_dofs = np.arange(1, self.nNodes - 1)

        # Material properties
        self.A_max = 0.01
        self.A_min = 0.001
        self.E = 70000e6                        # Young's Modulus

        # Dirichlet 
        self.U_m = 0.005                        # Amplitude of the prescribed displacement
        self.m = 2                              # Frequency
        self.time_span = 1                      # Time span
        self.steps = 10                         # Time steps
        self.vec_T = np.linspace(0, self.time_span, self.steps)
        self.U_t = self.U_m * np.sin(self.m * np.pi * self.vec_T / self.time_span)

        # Damage model
        self.sigma_u = 200e6                            # Ultimate stress
        self.r_0 = self.sigma_u/np.sqrt(self.E)         # Initial damage threshold
        self.beta_hs = 1.5
        self.q_inf = self.beta_hs*self.r_0
        self.A_hs = 0.1
        
        # Damage variables
        self.strain = np.zeros((self.nEle, self.steps))
        self.stress = np.zeros((self.nEle, self.steps))

        self.Etan = np.zeros((self.nEle, self.steps))

        self.r = np.zeros((self.nEle, self.steps))
        self.damage =np.zeros((self.nEle, self.steps))
  
class Solver_params:
    def __init__(self):
        self.max_iter = 10              # Number of maximum iteration by the solver
        self.tol = 1e-8                 # Minimum tolerance to achieve

class Reference_element:
    def __init__(self, p):
        if p == 1:
            self.X_ref = np.array([-1, 1])
            self.Z_g = 0
            self.nZ = 1
            self.W_g = 2

            self.Nxi = np.array([[-1/2, 1/2]])

def plot_displacement(Main_model):

    for i in range(Main_model.nNodes):
        plt.plot(Main_model.vec_T, Main_model.U[i, :], label=f'Node {i+1}')
    
    plt.xlabel('Position X (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacement at each element')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_strain(Main_model):
    for n in range(Main_model.steps):
        plt.plot(range(1, Main_model.nEle +1), Main_model.strain[:, n], label=f'Time Step {n+1}')
    
    plt.xlabel('Element')
    plt.ylabel('Strain')
    plt.title('Strain in each element over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_stress(Main_model):
    for n in range(Main_model.steps):
        plt.plot(range(1, Main_model.nEle +1), Main_model.stress[:, n], label=f'Time Step {n+1}')
    
    plt.xlabel('Element')
    plt.ylabel('Stress')
    plt.title('Stress in each element over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_strain_and_stress(Main_model):
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))  # Create two subplots, one for strain and one for stress

    # Plot strain
    for n in range(Main_model.steps):
        ax[0].plot(range(1, Main_model.nEle + 1), Main_model.strain[:, n], label=f'Time Step {n+1}')
    ax[0].set_xlabel('Element')
    ax[0].set_ylabel('Strain')
    ax[0].set_title('Strain over time')
    ax[0].legend()
    ax[0].grid(True)

    # Plot stress
    for n in range(Main_model.steps):
        ax[1].plot(range(1, Main_model.nEle + 1), Main_model.stress[:, n], label=f'Time Step {n+1}')
    ax[1].set_xlabel('Element')
    ax[1].set_ylabel('Stress (Pa)')
    ax[1].set_title('Stress over time')
    ax[1].legend()
    ax[1].grid(True)

    # Adjust the layout for better display
    plt.tight_layout()
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

    ani = FuncAnimation(fig, update, frames=Main_model.steps, interval=100)
    plt.show()

def Area(Main_model, x):
    A_max = Main_model.A_max
    A_min = Main_model.A_min

    L = Main_model.L

    if 0 <= x < 0.4*L or 0.6*L <= x <=L:
        A = A_max
    elif 0.4*L <= x <= 0.6*L:
        A = A_min

    return A

def q_hardening_law(Main_model,r):

    q_r = Main_model.q_inf -(Main_model.q_inf - Main_model.r_0)*np.exp(Main_model.A_hs*(1-r/Main_model.r_0))

    return q_r


def Compute_variables(Main_model, Element, U, n):
    # Compute strains and stresss of elements

    nEle = Main_model.nEle

    Nxi = Element.Nxi

    X = Main_model.X
    T = Main_model.T

    for ele in range(0,nEle):
        Te = T[ele,:]
        Xe = X[Te]
        h = Xe[1]- Xe[0]

        Nxi_ig = Nxi*2/h
        U_local = U[Te]

        d = Main_model.damage[ele,n]
        E = Main_model.E
        
        # Compute strain and stress
        strain =np.dot(Nxi_ig, U_local)
        stress = (1-d)*E*strain

        Main_model.strain[ele, n] = strain
        Main_model.stress[ele, n] = stress
    return print("Strains & Stress computed")


def Assemble_K_global(Main_model, Element):

    nNodes = Main_model.nNodes
    nEle = Main_model.nEle
    nZ = Element.nZ

    X = Main_model.X
    T = Main_model.T

    Nxi = Element.Nxi
    nZ = Element.nZ
    W_g = Element.W_g

    #Constituve model
    E = Main_model.E

    
    K_global = np.zeros((nNodes, nNodes))
    for ele in range(0,nEle):
        Te = T[ele,:]
        Xe = X[Te]

        
        X_midpoint = (Xe[1] + Xe[0])/2
        A = Area(Main_model, X_midpoint)

        h = Xe[1]- Xe[0]

        Ke = np.zeros((2,2))
        for ig in range(0, nZ):
            Nxi_ig = Nxi*2/h
            #print("Nxi", Nxi_ig)
            w_ig = W_g*h/2
            #print("Nxi_ig", Nxi_ig, "W_ig", w_ig)
            
            Ke += w_ig * (Nxi_ig.T @ (E * A * Nxi_ig))
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
        #print("Global K", K_global)

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
    #print("U_T", Main_model.U_t)

    U_step = np.zeros((Main_model.nNodes, 1), dtype=float)
    for n in range(Main_model.steps):
        print(f"Time step {n}")
        U_step[-1] = Main_model.U_t[n]
        #print("Imposed Ut", Main_model.U_t[n])

        # Solve nonlinear system
        U_step = Newton_Raphson(Main_model, Solver, Element, U_step, n)
        #print("U_step", U_step)

        Main_model.U[:,n] = U_step.flatten()  # Store results at each time step
        Compute_variables(Main_model, Element, Main_model.U[:,n], n)
        

    # Plot displacement, plot strains, plot stresses
    plot_displacement(Main_model)
    #plot_strain(Main_model)
    #plot_stress(Main_model)
    plot_strain_and_stress(Main_model)
    #animate_displacement(Main_model)
#-----------------------------------------------------------------


if __name__ == "__main__":
    main()