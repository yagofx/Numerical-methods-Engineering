import numpy as np
import csv
import matplotlib.pyplot as plt
import time as tp   
from matplotlib.animation import FuncAnimation

class Model_part:
    def __init__(self):

        #Time discretisation
        self.time_span = 1                                                  # Total time
        self.steps = 1000                                                   # Time steps
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


def main():

    start_time = tp.perf_counter()

    Main_model = Model_part()

    #Compute e11, e22, e33
    time_vector = Main_model.time_vector
    strain_tensor = Compute_strain(Main_model)
    plot_strain(time_vector, strain_tensor)


    end_time = tp.perf_counter()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")



#-----------------------------------------------------------------

if __name__ == "__main__":
    main()