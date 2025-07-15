#include <vector>
#include <array>

using Vec3 = std::array<double, 3>;

class TimeIntegrator {
public:
    std::vector<Vec3>& particle_coords;     // Reference to particle coords (updated in-place)
    const std::vector<Vec3>& velocity_points;  // Velocity sample points
    const std::vector<Vec3>& velocity_values;  // Velocity vectors at those points

    double dt;  // time step
    double T;   // total simulation time

    TimeIntegrator(std::vector<Vec3>& particle_coords_in,
                   const std::vector<Vec3>& vel_pts,
                   const std::vector<Vec3>& vel_vals,
                   double dt_in, double T_in);

    Vec3 InterpolationFunction(const Vec3& x) const;
    void ParticlesUpdate();

private:
    double RBF(const Vec3& x, const Vec3& xj) const;
};