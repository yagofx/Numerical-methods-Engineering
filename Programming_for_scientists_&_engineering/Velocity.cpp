#include "Velocity.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

void Velocity::ReadVelFile(const std::string& file)
{
    std::ifstream in(file);
    if(!in) throw std::runtime_error("Cannot open " + file);

    std::string line;
    std::size_t nPts = 0;

    /* ---------- 1) find and read the node coordinates ---------- */
    while(std::getline(in,line))
        if(line.rfind("POINTS",0)==0)          // “POINTS <N> float”
        {
            std::istringstream iss(line);
            std::string dummy; iss >> dummy >> nPts;
            break;
        }

    if(nPts==0) throw std::runtime_error("POINTS section not found in "+file);

    X_Gauss.resize(nPts);
    for(std::size_t i=0;i<nPts;++i)
        in >> X_Gauss[i][0] >> X_Gauss[i][1] >> X_Gauss[i][2];

    /* ---------- 2) skip to the velocity field  ----------------- */
    // For UNSTRUCTURED_GRID we want
    //   POINT_DATA <N>
    //   VECTORS   <name> <float/double>
    while(std::getline(in,line))
        if(line.rfind("VECTORS",0)==0) break;

    if(in.fail())
        throw std::runtime_error("VECTORS field not found in "+file);

    /* ---------- 3) read the vectors ---------------------------- */
    vels_Gauss.resize(nPts);
    for(std::size_t i=0;i<nPts;++i)
        in >> vels_Gauss[i][0] >> vels_Gauss[i][1] >> vels_Gauss[i][2];

    if(in.fail())
        throw std::runtime_error("Error while reading velocity vectors in "+file);
}