#include <vector>
#include <array>
#include <string>

using Vec3 = std::array<double, 3>;

class Velocity {
public:
    std::vector<Vec3> X_Gauss;
    std::vector<Vec3> vels_Gauss;

    void ReadVelFile(const std::string& file);
};