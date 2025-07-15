#include "Particles.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

void Particles::ReadParticlesFile(const std::string& file) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("Cannot open " + file);

    std::string line;
    size_t nPts = 0;

    while (std::getline(in, line)) {
        if (line.rfind("POINTS", 0) == 0) {
            std::istringstream iss(line);
            std::string dummy; iss >> dummy >> nPts;
            break;
        }
    }

    coords.resize(nPts);
    for (size_t i = 0; i < nPts; ++i)
        in >> coords[i][0] >> coords[i][1] >> coords[i][2];

    while (std::getline(in, line)) {
        if (line.rfind("SCALARS", 0) == 0) break;
    }
    std::getline(in, line); // skip LOOKUP_TABLE

    ids.resize(nPts);
    for (size_t i = 0; i < nPts; ++i)
        in >> ids[i];
}

void Particles::WriteParticlesFile(const std::string& file) const {
    std::ofstream out(file);
    if (!out) throw std::runtime_error("Cannot write " + file);

    out << "# vtk DataFile Version 3.0\nParticles\nASCII\n\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << coords.size() << " float\n";
    for (const auto& p : coords)
        out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';

    out << "\nVERTICES " << coords.size() << ' ' << coords.size() * 2 << '\n';
    for (size_t i = 0; i < coords.size(); ++i)
        out << "1 " << i << '\n';

    out << "\nPOINT_DATA " << coords.size() << '\n';
    out << "SCALARS myParticles int\nLOOKUP_TABLE T1\n";
    for (int id : ids)
        out << id << '\n';
}