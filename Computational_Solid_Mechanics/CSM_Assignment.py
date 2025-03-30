#-----------------------------------------------------------------
#
#               Computational Solid Mechanics - Assignment damage
#                   
#               Authors: Gabriel Ayu, David Conteras, Yago Trias
#   
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#
#               FE for 1D Beam damage model
#
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import time as tp   
from matplotlib.animation import FuncAnimation

#-----------------------------------------------------------------
#                       Classes
#-----------------------------------------------------------------

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
        self.A_max = 0.01
        self.A_min = 0.001
        self.E = 70000e6                        # Young's Modulus

        # Dirichlet 
        self.U_m = 0.005                        # Amplitude of the prescribed displacement
        self.m = 2                              # Frequency
        self.time_span = 1                      # Time span
        self.steps = 40                         # Time steps
        self.vec_T = np.linspace(0, self.time_span, self.steps)
        self.U_t = self.U_m * np.sin(self.m * np.pi * self.vec_T / self.time_span)

        # Damage model
        self.sigma_u = 200e6                            # Ultimate stress
        self.r_0 = self.sigma_u/np.sqrt(self.E)         # Initial damage threshold
        self.beta_hs = 1.5
        self.q_inf = self.beta_hs*self.r_0
        self.A_hs = 0.1

        # Damage variables
        self.r = np.zeros((self.nEle, self.steps))
        self.r[:,0] = self.r_0
        self.damage =np.zeros((self.nEle, self.steps))
        
        self.strain = np.zeros((self.nEle, self.steps))
        self.stress = np.zeros((self.nEle, self.steps))

        self.Etan = np.zeros((self.nEle, self.steps+1))
        self.Etan[:,0] = self.E

  
class Solver_params:
    def __init__(self):
        self.max_iter = 10              # Number of maximum iteration by the solver
        self.tol = 1e-16                # Minimum tolerance to achieve

class Reference_element:
    def __init__(self, p):
        if p == 1:
            self.X_ref = np.array([-1, 1])
            self.Z_g = 0
            self.nZ = 1
            self.W_g = 2

            self.Ni = np.array([[(1-self.Z_g)/2, (1+self.Z_g)/2,]])
            self.Nxi = np.array([[-1/2, 1/2]])


#-----------------------------------------------------------------
#                   Plotting
#-----------------------------------------------------------------
def plot_displacement(Main_model):
    for n in range(Main_model.steps):
        plt.plot(Main_model.X, Main_model.U[:, n], label=f'Time Step {n+1}')
    
    plt.xlabel('Position X (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacement at each node over time')
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
    for n in range(0,Main_model.steps,1):
        ax[0].plot(range(1, Main_model.nEle + 1), Main_model.strain[:, n], label=f'Time Step {n+1}')
    ax[0].set_xlabel('Element')
    ax[0].set_ylabel('Strain')
    ax[0].set_title('Strain over time')
    ax[0].legend()
    ax[0].grid(True)

    # Plot stress
    for n in range(0,Main_model.steps,1):
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

def plot_midpoint_data(Main_model):
    midpoint_element = 2  # Adjust for zero-based indexing

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

    # Plot stress over time
    axes[0].plot(Main_model.vec_T, Main_model.stress[midpoint_element, :], label='Stress')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('Stress')
    axes[0].set_title(f'Stress at Element {midpoint_element + 1}')
    axes[0].legend()
    axes[0].grid(True)

    # Plot damage over time
    axes[1].plot(Main_model.vec_T, Main_model.damage[midpoint_element, :], label='Damage', color='red')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Damage')
    axes[1].set_title(f'Damage at Element {midpoint_element + 1}')
    axes[1].legend()
    axes[1].grid(True)

    # Plot stress vs strain
    axes[2].plot(Main_model.strain[midpoint_element, :], Main_model.stress[midpoint_element, :], label='Stress-Strain', color='orange')
    axes[2].set_xlabel('Strain')
    axes[2].set_ylabel('Stress')
    axes[2].set_title(f'Stress-Strain at Element {midpoint_element + 1}')
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()


#-----------------------------------------------------------------
#                   FE Model
#-----------------------------------------------------------------

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

    q_inf = Main_model.q_inf
    r_0 = Main_model.r_0
    A_hs = Main_model.A_hs

    q_r = q_inf -(q_inf - r_0)*np.exp(A_hs*(1-r/r_0))

    return q_r

def H_hardening_law(Main_model,r):

    q_inf = Main_model.q_inf
    r_0 = Main_model.r_0
    A_hs = Main_model.A_hs

    H_r = A_hs*(q_inf - r_0)/r_0*np.exp(A_hs*(1-r/r_0))

    return H_r


def Update_damage(Main_model, n):
    # Update damage model

    nEle = Main_model.nEle
    E = Main_model.E

    for ele in range(0,nEle):
        strain = Main_model.strain[ele, n]
        if n == 0:
            r = Main_model.r_0
            damage = 0
        else:
            r = Main_model.r[ele,n-1]
            damage = Main_model.damage[ele, n-1]

        # Step 1
        sigma_bar = E*strain
        tau_epsilon = np.sqrt(strain*sigma_bar)


        # Step 2
        if tau_epsilon <= r:
            r_updated = r
            damage_updated = damage
            stress_updated = (1-damage_updated)*sigma_bar
            Etan_updated = (1-damage_updated)*E 


        # Step 3
        elif tau_epsilon > r:
            r_updated = tau_epsilon
            q = q_hardening_law(Main_model, r_updated)
            H = H_hardening_law(Main_model,r_updated)
            damage_updated = 1-q/r_updated
            stress_updated = (1-damage_updated)*sigma_bar
            Etan_updated = (1-damage_updated)*E - ((q-H*r_updated)/(r_updated**3))*(sigma_bar*sigma_bar)

        #Update values
        #print("strain:", strain,"sigma_bar", sigma_bar,"tau_epsi", tau_epsilon, "r", r_updated, "damage: ", damage_updated, )
        Main_model.r[ele, n] = r_updated
        Main_model.damage[ele, n] = damage_updated
        Main_model.stress[ele, n] = stress_updated
        Main_model.Etan[ele, n+1] = Etan_updated

    return print("Damage model updated")
    

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
        U_local = U[Te].reshape(-1, 1)

        # Compute strain
        strain =np.dot(Nxi_ig, U_local)

        Main_model.strain[ele, n] = strain
    return print("Strains computed")


def Assemble_K_global(Main_model, Element, n):

    nNodes = Main_model.nNodes
    nEle = Main_model.nEle
    nZ = Element.nZ

    X = Main_model.X
    T = Main_model.T

    Ni = Element.Ni
    Nxi = Element.Nxi
    nZ = Element.nZ
    W_g = Element.W_g
    
    K_global = np.zeros((nNodes, nNodes))
    for ele in range(0,nEle):
        Te = T[ele,:]
        Xe = X[Te]

        Etan = Main_model.Etan[ele,n]

        h = Xe[1]- Xe[0]

        Ke = np.zeros((2,2))
        for ig in range(0, nZ):
            Ni_ig = Ni
            Nxi_ig = Nxi*2/h
            w_ig = W_g*h/2

            X_ig = np.dot(Ni_ig, Xe)
            A_ig = Area(Main_model, X_ig)
            Etan = Main_model.Etan[ele,n]

            Ke += A_ig*w_ig * (Nxi_ig.T @ (Etan * Nxi_ig))

        K_global[np.ix_(Te, Te)] += Ke


    return K_global

def Newton_Raphson(Main_model, Solver, Element, U_step, n):

    max_iter = Solver.max_iter
    tol = Solver.tol

    fixed_dofs = Main_model.fixed_dofs
    free_dofs = Main_model.free_dofs

    K_global = Assemble_K_global(Main_model, Element, n)

    for i in range(max_iter + 1):

        R = -np.dot(K_global, U_step)
        Rnorm = np.linalg.norm(R[free_dofs])
        print(f"NL iter: {i}, |R|: {Rnorm}")

        if Rnorm < tol:
            print(f"Converged in {i} iterations!")
            return U_step
        
    
        K_uu = K_global[free_dofs, :][:, free_dofs]
        
        delta_u = np.linalg.solve(K_uu, R[free_dofs])
        
        U_step[free_dofs] += delta_u
        U_step[fixed_dofs] = np.array([0.0, Main_model.U_t[n]]).reshape(-1, 1)
        #print(f"U {i}", U_step)

    return U_step

def main():

    start_time = tp.perf_counter()

    Main_model = Model_part()
    Solver = Solver_params()

    p = 1
    Element = Reference_element(p)

    Main_model.U = np.zeros((Main_model.nNodes, Main_model.steps), dtype=float)

    U_step = np.zeros((Main_model.nNodes, 1), dtype=float)
    for n in range(Main_model.steps):
        print("%-----------------------------------------------------------------------")
        print(f"Time step {n}")
        print("%-----------------------------------------------------------------------")
        U_step[-1] = Main_model.U_t[n]


        # Solve nonlinear system

        U_step = Newton_Raphson(Main_model, Solver, Element, U_step, n)

        Main_model.U[:,n] = U_step.flatten()
        Compute_variables(Main_model, Element, Main_model.U[:,n], n)
        Update_damage(Main_model, n)
        

    end_time = tp.perf_counter()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")

    plot_displacement(Main_model)
    plot_strain_and_stress(Main_model)
    plot_midpoint_data(Main_model)
    #animate_displacement(Main_model)
#-----------------------------------------------------------------


if __name__ == "__main__":
    main()