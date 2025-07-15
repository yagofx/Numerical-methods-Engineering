#include "TimeIntegrator.h"
#include <cmath>

TimeIntegrator::TimeIntegrator(std::vector<Vec3>& particle_coords_in,
                               const std::vector<Vec3>& vel_pts,
                               const std::vector<Vec3>& vel_vals,
                               double dt_in, double T_in)
    : particle_coords(particle_coords_in), velocity_points(vel_pts),
      velocity_values(vel_vals), dt(dt_in), T(T_in) {}

double TimeIntegrator::RBF(const Vec3& x, const Vec3& xj) const {
    const double epsilon = 10;    // your fixed epsilon
    const double r_cut = 0.5;      // your chosen cutoff radius

    double dist2 = 0.0;
    for (int i = 0; i < 3; ++i)
        dist2 += (x[i] - xj[i]) * (x[i] - xj[i]);

    double dist = std::sqrt(dist2);

    if (dist > r_cut)
        return 0.0;  // zero weight outside cutoff radius

    return std::exp(-epsilon * epsilon * dist2);
}

Vec3 TimeIntegrator::InterpolationFunction(const Vec3& x) const {
    Vec3 interpolated_vel = {0.0, 0.0, 0.0};
    double w_sum = 0.0;

    for (size_t j = 0; j < velocity_points.size(); ++j) {
        double w = RBF(x, velocity_points[j]);
        for (int i = 0; i < 3; ++i)
            interpolated_vel[i] += w * velocity_values[j][i];
        w_sum += w;
    }

    if (w_sum > 1e-12) {
        for (int i = 0; i < 3; ++i)
            interpolated_vel[i] /= w_sum;
    }

    return interpolated_vel;
}

void TimeIntegrator::ParticlesUpdate() {
    for (auto& x : particle_coords) {
        Vec3 v = InterpolationFunction(x);
        for (int i = 0; i < 3; ++i)
            x[i] += dt * v[i];
    }
}