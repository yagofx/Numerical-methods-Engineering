#include <iostream>
#include <iomanip>
#include <sstream>
#include "Velocity.h"
#include "Particles.h"
#include "TimeIntegrator.h"

int main()
{
    try {
        // Load velocity field
        Velocity velocity;
        velocity.ReadVelFile("velocity.vtk");

        // Load initial particles
        Particles particles;
        particles.ReadParticlesFile("particles.0001.vtk");

        // Time settings
        double T = 4;        // total simulation time
        double delta_t = 0.001; // time step
        int num_steps = static_cast<int>(T / delta_t);

        // Create TimeIntegrator with references to velocity points/values and particle coordinates
        TimeIntegrator integrator(particles.coords, velocity.X_Gauss, velocity.vels_Gauss, delta_t, T);

        for (int step = 0; step <= num_steps; ++step) {
            if (step == num_steps)
                break; // finished simulation

            // Update particle positions for next timestep
            integrator.ParticlesUpdate();

            std::ostringstream filename;
            filename << "output/particles_t" << std::setfill('0') << std::setw(4) << step << ".vtk";
            particles.WriteParticlesFile(filename.str());
        }

    } catch(const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}