#-----------------------------------------------------------------
#
#               Computational Solid Mechanics - Assignment damage
#                   
#               Authors: Gabriel Ayu, Yago Trias
#   
#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#
#               FE for 1D Beam damage model
#
#-----------------------------------------------------------------

import numpy as np
import csv
import matplotlib.pyplot as plt
import time as tp   
from matplotlib.animation import FuncAnimation
import generate_VTK

#-----------------------------------------------------------------
#                       Classes
#-----------------------------------------------------------------

class Model_part:
    def __init__(self):
        self.L = 1                                                  # Length of the beam (m)
        self.nEle = 5                                              # Number of elements
        self.nNodes = self.nEle + 1                                 # Number of nodes
        self.X = np.linspace(0,self.L,self.nNodes)                  # Coordinates of nodes 
        self.T = np.array([[i, i + 1] for i in range(self.nEle)])   # Connectivity of the elements

        self.fixed_dofs = np.array([0,self.nNodes-1])
        self.free_dofs = np.arange(1, self.nNodes - 1)

        # Material properties
        self.A_max = 0.01
        self.A_min = 0.001

        # Dirichlet 
        self.U_m = 0.001                                # Amplitude of the prescribed displacement
        self.m = 2.5                                    # Frequency
        self.time_span = 5                              # Time span
        self.steps = 500 + 1                            # Time steps
        self.delta_t = self.time_span/self.steps        # Delta t
        self.vec_T = np.linspace(0, self.time_span, self.steps)

        self.U_t = self.U_m * np.sin(self.m * np.pi * self.vec_T / self.time_span)
        #self.U_t = self.U_m * self.vec_T * np.sin(self.m * np.pi * self.vec_T / self.time_span)

        # Plasticity model
        self.E = 2.1e+11            # Young's Modulus
        self.K = 5.0e+10            # Isotropic hardening parameters
        #self.H = 1.0e+10            # Kinematic hardening parameters
        self.eta = 5.0e+10          # Viscosity parameter
        # self.K = 0
        self.H = 0
        self.sigma_y = 4.2e+8       # Yield Stress
        self.sigma_inf = 2.0e+9     # Maximum Asymptotic Stress
        self.delta = 150.0          # Saturation parameter


        self.C = np.diag([self.E, self.K, self.H])

        self.Etan = self.E*np.ones((self.nEle,1))   # Ctang

        self.strain = np.zeros((self.nEle))
        self.stress = np.zeros((self.nEle,3))

        # Plasticity variables
        self.E_p = np.zeros((self.nEle,3))      # Plastic strain

        #Stored stress & strain
        self.strain2 = np.zeros((self.nEle, self.steps))
        self.stress2 = np.zeros((self.nEle, self.steps))
        self.plasticity = np.zeros((self.nEle, self.steps))



  
class Solver_params:
    def __init__(self):
        self.max_iter = 10             # Number of maximum iteration by the solver
        self.tol = 1e-8                # Minimum tolerance to achieve

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
    for n in range(0,Main_model.steps,20):
        plt.plot(Main_model.X, Main_model.U[:, n], label=f'Time Step {n}')
    
    plt.xlabel('Position X (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacement at each node over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_strain(Main_model):
    for n in range(Main_model.steps):
        plt.plot(range(1, Main_model.nEle +1), Main_model.strain2[:, n], label=f'Time Step {n}')
    
    plt.xlabel('Element')
    plt.ylabel('Strain')
    plt.title('Strain in each element over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_stress(Main_model):
    for n in range(Main_model.steps):
        plt.plot(range(1, Main_model.nEle +1), Main_model.stress2[:, n], label=f'Time Step {n}')
    
    plt.xlabel('Element')
    plt.ylabel('Stress')
    plt.title('Stress in each element over time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_strain_and_stress(Main_model):
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))  # Create two subplots, one for strain and one for stress

    # Plot strain
    for n in range(0,Main_model.steps,20):
        ax[0].plot(range(1, Main_model.nEle + 1), Main_model.strain2[:, n], label=f'Time Step {n}')
    ax[0].set_xlabel('Element')
    ax[0].set_ylabel('Strain')
    ax[0].set_title('Strain over time')
    ax[0].legend()
    ax[0].grid(True)

    # Plot stress
    for n in range(0,Main_model.steps,20):
        ax[1].plot(range(1, Main_model.nEle + 1), Main_model.stress2[:, n], label=f'Time Step {n}')
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
    midpoint_element = round(Main_model.nEle/2)  # Adjust for zero-based indexing

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

    # Plot stress over time
    axes[0].plot(Main_model.vec_T, Main_model.stress2[midpoint_element, :], label='Stress')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('Stress')
    axes[0].set_title(f'Stress at Element {midpoint_element + 1}')
    axes[0].legend()
    axes[0].grid(True)

    # Plot plasticity over time
    axes[1].plot(Main_model.vec_T, Main_model.plasticity[midpoint_element, :], label='Plasticity', color='red')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Plasticity')
    axes[1].set_title(f'Plasticity at Element {midpoint_element + 1}')
    axes[1].legend()
    axes[1].grid(True)

    # Plot stress vs strain
    axes[2].plot(Main_model.strain2[midpoint_element, :], Main_model.stress2[midpoint_element, :], label='Stress-Strain', color='orange')
    axes[2].set_xlabel('Strain')
    axes[2].set_ylabel('Stress')
    axes[2].set_title(f'Stress-Strain at Element {midpoint_element + 1}')
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()

def export_results_to_csv(Main_model, prefix="results"):
    """
    Exports r, damage, strain, and stress matrices to CSV files.
    Each file will have shape (nEle x (steps + 1)), with rows = elements and columns = [x_center, step_0, ..., step_N].
    """
    variables = {
        "r": Main_model.r,
        "damage": Main_model.damage,
        "strain": Main_model.strain,
        "stress": Main_model.stress
    }

    # Compute element x-centers (mean of node coordinates)
    x_coords = Main_model.X.flatten()
    x_centers = np.mean(x_coords[Main_model.T], axis=1).reshape(-1, 1)  # (nEle, 1)

    for name, matrix in variables.items():
        data = np.hstack((x_centers, matrix))  # shape (nEle, steps+1)
        header = "x," + ",".join([f"step_{i}" for i in range(Main_model.steps)])
        filename = f"{prefix}_{name}.csv"
        np.savetxt(filename, data, delimiter=",", header=header, comments='')
        print(f"Saved {name} to {filename}")
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

def pi_prime(Main_model, x):

    sigma_inf = Main_model.sigma_inf
    sigma_y = Main_model.sigma_y
    delta = Main_model.delta
    K = Main_model.K

    output = (sigma_inf - sigma_y)*(1 - np.exp(-delta*x)) + K*x

    return output

def pi_pprime(Main_model, x):
    
    sigma_inf = Main_model.sigma_inf
    sigma_y = Main_model.sigma_y
    delta = Main_model.delta

    output = delta*(sigma_inf - sigma_y)*np.exp(-delta*(x))

    return output


def Update_Plasticity(Main_model, Solver, Element, U_step):
    nEle = Main_model.nEle
    E = Main_model.E
    K = Main_model.K
    H = Main_model.H
    sigma_y = Main_model.sigma_y
    sigma_inf = Main_model.sigma_inf
    delta = Main_model.delta

    delta_t = Main_model.delta_t
    eta = Main_model.eta

    E_p = Main_model.E_p
    strain = Compute_Strain(Main_model, Element, U_step)
    stress = Main_model.stress
    Main_model.strain = strain
    #print("Total strain:", strain)

    for ele in range(nEle):

        sstress, qq, qq_bar = stress[ele]

        strain_local = strain[ele]
        strain_p, chi, chi_bar = E_p[ele]

        strain_e = strain_local - strain_p
        sigma = strain_e * E
        q = -pi_prime(Main_model, chi)
        q_bar = -chi_bar * H

        f_trial = np.abs(sigma - q_bar) - sigma_y + q

        if f_trial <= 0:
            Etan_updated = E

            sigma_updated = sigma
            q_updated = qq
            q_bar_updated = qq_bar

            strain_p_updated = strain_p
            chi_updated = chi
            chi_bar_updated = chi_bar

        else:
            gamma = Newton_Raphson_gamma(Main_model, Solver, f_trial, chi)
            plastic_multi = gamma*delta_t
            sign = np.sign(sigma - q_bar)

            sigma_updated = sigma - plastic_multi*E*sign
            pi_p = pi_prime(Main_model, chi + gamma*delta_t)  
            q_updated = -pi_p
            q_bar_updated = q_bar + plastic_multi*H*sign
            
            pi_pp = pi_pprime(Main_model, chi + gamma*delta_t) 
            Etan_updated = E * (1 - E*(E + pi_pp + H + eta/delta_t)**(-1))

            strain_p_updated = strain_p + plastic_multi * sign
            chi_updated = chi + plastic_multi
            chi_bar_updated = chi_bar - plastic_multi * sign

        stress_updated = np.array([[sigma_updated, q_updated, q_bar_updated]])
        strain_p_vec = np.array([[strain_p_updated, chi_updated, chi_bar_updated]])

        print(f"ele: {ele + 1}, strain_e: {strain_e}, f_trial: {f_trial}, sigma: {sigma}, sign: {np.sign(sigma - q_bar)}, gamma: {gamma if f_trial > 0 else 0}, strain_p_upd: {strain_p_updated}")

        # Save updated variables
        Main_model.E_p[ele] = strain_p_vec
        Main_model.stress[ele, :] = stress_updated
        Main_model.Etan[ele] = Etan_updated

    return

def Compute_Strain(Main_model, Element, U):
    # Compute strains and stresss of elements

    nEle = Main_model.nEle

    Nxi = Element.Nxi

    X = Main_model.X
    T = Main_model.T

    strain = np.zeros(nEle)
    for ele in range(nEle):
        Te = T[ele,:]
        Xe = X[Te]
        h = Xe[1]- Xe[0]

        Nxi_ig = Nxi*2/h
        U_local = U[Te].reshape(-1, 1)

        # Compute strain
        strain_local = np.dot(Nxi_ig, U_local)

        strain[ele] = strain_local
        
    return strain

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
    
    K_global = np.zeros((nNodes, nNodes))
    for ele in range(0,nEle):
        Te = T[ele,:]
        Xe = X[Te]

        Etan = Main_model.Etan[ele]

        h = Xe[1]- Xe[0]

        Ke = np.zeros((2,2))
        for ig in range(0, nZ):
            Ni_ig = Ni
            Nxi_ig = Nxi*2/h
            w_ig = W_g*h/2

            X_ig = np.dot(Ni_ig, Xe)
            A_ig = Area(Main_model, X_ig)

            Ke += A_ig*w_ig * (Nxi_ig.T @ (Etan * Nxi_ig))

        K_global[np.ix_(Te, Te)] += Ke
    return K_global

def f_internal(Main_model, Element, U_step):
    nNodes = Main_model.nNodes
    nEle = Main_model.nEle
    nZ = Element.nZ

    X = Main_model.X
    T = Main_model.T

    Ni = Element.Ni
    Nxi = Element.Nxi
    W_g = Element.W_g

    sigma = Main_model.stress[:,0].copy()

    f_int = np.zeros((nNodes, 1))

    for ele in range(nEle):
        Te = T[ele, :]
        Xe = X[Te]
        h = Xe[1] - Xe[0]

        stress_local = sigma[ele]

        f_ele = np.zeros((2, 1))
        for ig in range(nZ):
            Ni_ig = Ni
            Nxi_ig = Nxi[ig, :] * 2 / h
            Nxi_ig = Nxi_ig[:, None]
            w_ig = W_g * h / 2 

            X_ig = np.dot(Ni_ig, Xe)
            A_ig = Area(Main_model, X_ig)

            sigma_ig = stress_local

            f_ele += A_ig * w_ig * sigma_ig * Nxi_ig

        f_int[Te] += f_ele

    return f_int

def Newton_Raphson(Main_model, Solver, Element, U_step, n):

    max_iter = Solver.max_iter
    tol = Solver.tol

    fixed_dofs = Main_model.fixed_dofs
    free_dofs = Main_model.free_dofs

    for i in range(max_iter + 1):
        R = -f_internal(Main_model, Element, U_step)
        Rnorm = np.linalg.norm(R[free_dofs], "fro")
        print(f"\nNL iter: {i + 1}, |R|: {Rnorm}")

        if Rnorm < tol and i != 0:
            print(f"\nConverged in {i + 1} iterations!\n")
            return U_step

        K_global = Assemble_K_global(Main_model, Element)
        K_uu = K_global[free_dofs, :][:, free_dofs]

        delta_u = np.linalg.solve(K_uu, R[free_dofs])
        U_step[free_dofs] += delta_u
        U_step[fixed_dofs] = np.array([0.0, Main_model.U_t[n]]).reshape(-1, 1)

        Update_Plasticity(Main_model, Solver ,Element, U_step)

    return U_step

def Newton_Raphson_gamma(Main_model, Solver, f_trial, chi):

    max_iter = Solver.max_iter
    tol = Solver.tol

    E = Main_model.E
    K = Main_model.K
    H = Main_model.H
    eta = Main_model.eta
    delta_t = Main_model.delta_t

    gamma = 0

    for i in range(max_iter + 1):
        pi_p_chi_gamma = pi_prime(Main_model, chi + gamma*delta_t)
        pi_p_chi = pi_prime(Main_model, chi)            
        g = f_trial - gamma*delta_t*(E + H + eta/delta_t) - (pi_p_chi_gamma - pi_p_chi)

        pi_pp = pi_pprime(Main_model, chi + gamma*delta_t)
        K = -(E + pi_pp + H + eta/delta_t)*delta_t                                                                # Derivative of g

        gamma_updated = gamma - g/K
        gamma = gamma_updated
        R = np.abs(g)

        print(f" gamma NL iter: {i + 1}, |R|: {R}")

        if R < tol:
            print(f"\n Converged in {i + 1} iterations!\n")
            return gamma

    if R > tol:
        print("\n    Warning no convergence.\n")

    return gamma

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
        print("U_step", U_step)


        # Solve nonlinear system
        U_step = Newton_Raphson(Main_model, Solver, Element, U_step, n)
        Main_model.U[:,n] = U_step.flatten()
        Main_model.strain2[:, n] = Main_model.strain[:]
        Main_model.stress2[:, n] = Main_model.stress[:,0]
        Main_model.plasticity[:, n] = Main_model.E_p[:,0]
        #generate_VTK.write_vtk(Main_model, n, filename="solution")
        

    end_time = tp.perf_counter()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")

    # plot_displacement(Main_model)
    # plot_strain_and_stress(Main_model)
    plot_midpoint_data(Main_model)
    #export_results_to_csv(Main_model)

    data = np.vstack((Main_model.strain2[2,:], Main_model.stress2[2,:])).T  # shape: (steps, 2)
    np.savetxt("NL_hardening_RateD.csv", data, delimiter=',', header='strain,stress', comments='')

    #animate_displacement(Main_model)
#-----------------------------------------------------------------


if __name__ == "__main__":
    main()