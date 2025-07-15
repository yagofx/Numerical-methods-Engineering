#include <vector>
#include <array>
#include <string>

using Vec3 = std::array<double, 3>;

class Particles {
public:
    std::vector<Vec3> coords;
    std::vector<int>  ids;

    void ReadParticlesFile(const std::string& file);
    void WriteParticlesFile(const std::string& file) const;
};