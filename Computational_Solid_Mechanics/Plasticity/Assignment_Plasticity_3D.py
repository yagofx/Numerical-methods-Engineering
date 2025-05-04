import numpy as np
import csv
import matplotlib.pyplot as plt
import time as tp   
from matplotlib.animation import FuncAnimation

class Model_part:
    def __init__(self):

        #Time discretisation
        self.time_span = 1                                                  # Total time
        self.steps = 200                                                   # Time steps
        self.delta_t = self.time_span/self.steps                            # delta t
        self.time_vector = np.linspace(0, self.time_span, self.steps)       # time vector

        #Input data
        self.input_e11 = np.array([[0, 0.02, -0.02, 0.02],[0, self.time_span/3, 2*self.time_span/3, self.time_span]])             # Loading path
           

        # Plasticity model
        self.E = 2.1e+11            # Young's Modulus
        self.sigma_y = 4.2e+8       # Yield Stress
        self.sigma_inf = 2.0e+9     # Maximum asymptotic stress

        self.K = 5.0e+10            # Isotropic hardening parameter
        self.H = 1.0e+10            # Kinematic hardening parameter
        self.eta = 5.0e+10          # Viscosity parameter
        self.nu = 0.3               # Poisson's coefficient
        self.delta = 150            # Saturation parameter

        # Lame parameters
        self.lamb = self.nu*self.E/((1 + self.nu)*(1 - 2*self.nu))  # Lame lambda
        self.mu = self.E/(2*(1 + self.nu))  # Shear elastic modulus

        self.C = np.zeros((3, 3, 3, 3))

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        delta_ij = int(i == j)
                        delta_kl = int(k == l)
                        delta_ik = int(i == k)
                        delta_jl = int(j == l)
                        delta_il = int(i == l)
                        delta_jk = int(j == k)

                        self.C[i, j, k, l] = (
                            self.lamb * delta_ij * delta_kl +
                            self.mu * (delta_ik * delta_jl + delta_il * delta_jk)
                        )

        #Internal variables
        self.strain_p = np.zeros((3,3))     #3x3
        self.chi = 0.0                      #scalar 
        self.chi_bar = np.zeros((3,3))      #3x3

        self.stress = np.zeros((3,3))       #stress tensor
        self.q = 0.0                        #scalar
        self.q_bar = np.zeros((3,3))        #3x3       

        #Stored stress & strain

def Compute_strain(Main_model):

    time_vector = Main_model.time_vector

    input_strain = Main_model.input_e11[0, :]
    input_time = Main_model.input_e11[1, :]

    strain_e11 = np.interp(time_vector, input_time, input_strain)

    nu = Main_model.nu
    strain_e22 = -nu*strain_e11
    strain_e33 = -nu*strain_e11

    strain_tensor = np.vstack((strain_e11, strain_e22, strain_e33)).T


    return strain_tensor

def plot_strain(time, strain):
    plt.plot(time, strain[:, 0], label='e11')
    plt.plot(time, strain[:, 1], label='e22')
    plt.plot(time, strain[:, 2], label='e33')
    plt.xlabel('Time (s)')
    plt.ylabel('Strain')
    plt.title('Strain components over time')
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_stress_strain(strain11, stress11):
    plt.plot(strain11, stress11, label=r'$\sigma_{11}$ vs $\varepsilon_{11}$')
    plt.xlabel(r'Strain $\varepsilon_{11}$')
    plt.ylabel(r'Stress $\sigma_{11}$ (Pa)')
    plt.title('Uniaxial Stress-Strain Response')
    plt.grid(True)
    plt.legend()
    plt.show()

def Update_Plasticity(Main_model, input_strain_tensor):

    C = Main_model.C
    K = Main_model.K
    H = Main_model.H
    sigma_y = Main_model.sigma_y
    mu = Main_model.mu

    strain_p = Main_model.strain_p
    chi = Main_model.chi
    chi_bar = Main_model.chi_bar

    total_strain = input_strain_tensor
    strain_e = total_strain - strain_p

    sigma_trial = np.einsum('ijkl,kl->ij', C, strain_e)
    q_trial = -K * chi
    q_bar_trial = - (2.0 / 3.0)*H*chi_bar 

    # Deviatoric part of trial stress
    sigma_m = np.trace(sigma_trial)/3.0
    sigma_dev = sigma_trial - sigma_m * np.eye(3)

    # Yield function
    s_rel = sigma_dev - q_bar_trial
    norm_s = np.linalg.norm(s_rel)
    f_trial = norm_s - np.sqrt(2.0 / 3.0)*(sigma_y - q_trial)

    print("f_trial:", f_trial)

    if f_trial <= 0:

        gamma = 0.0
        n = np.zeros((3,3))

        sigma_updated = sigma_trial
        q_updated = q_trial
        q_bar_updated = q_bar_trial

        strain_p_updated = strain_p
        chi_updated = chi
        chi_bar_updated = chi_bar

    else:

        gamma = (2*mu + 2/3*K + 2/3*H)**(-1)*f_trial

        n = s_rel / norm_s

        sigma_updated = sigma_trial - 2*mu*gamma*n
        q_updated = q_trial - K*gamma*np.sqrt(2.0 / 3.0)
        q_bar_updated = q_bar_trial + (2.0 / 3.0)*H*gamma*n

        strain_p_updated = strain_p + gamma*n
        chi_updated = chi + gamma*np.sqrt(2.0 / 3.0)
        chi_bar_updated = chi_bar - gamma*n
    
    #Update variables in object
    Main_model.stress = sigma_updated
    Main_model.q = q_updated
    Main_model.q_bar = q_bar_updated

    Main_model.strain_p = strain_p_updated
    Main_model.chi = chi_updated
    Main_model.chi_bar = chi_bar_updated

def main():

    start_time = tp.perf_counter()

    Main_model = Model_part()

    #Compute e11, e22, e33
    time_vector = Main_model.time_vector
    strain = Compute_strain(Main_model)
    input_strain = np.column_stack((time_vector, strain))  # shape: (steps, 4)
    np.savetxt("input_strain.csv", input_strain, delimiter=",", 
           header="time,e11,e22,e33", comments="")
    plot_strain(time_vector, strain)

    stress = np.zeros((Main_model.steps))
    for n in range(Main_model.steps):
        print("%-----------------------------------------------------------------------")
        print(f"Time step {n}")
        print("%-----------------------------------------------------------------------")

        #Compute strain tensor
        strain_vector = strain[n,:]
        print("strain tensor", strain_vector)

        #Update plasticity function
        strain_tensor = np.diag(strain_vector)
        Update_Plasticity(Main_model, strain_tensor)

        stress[n] = Main_model.stress[0,0]  #save stress11
        


        



    end_time = tp.perf_counter()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")

    plot_stress_strain(strain[:, 0], stress)



#-----------------------------------------------------------------

if __name__ == "__main__":
    main()