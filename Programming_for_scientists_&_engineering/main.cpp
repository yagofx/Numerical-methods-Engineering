#include <iostream>
#include <fstream>
#include <string>
using namespace std;

const int nNODES = 3917;
const int nELEM = 7182;

int main()
{   
    //---------------------------------------
    //      Read mesh
    //---------------------------------------
    const int col = 3;
    double X[nNODES][col];  // Coordinates
    int T[nELEM][col];      // Connectivity
    double U[nNODES];       // Nodal solution (scalar)

    ifstream ob("mesh.vtk");
    if (!ob) {
        cerr << "Cannot open file mesh.vtk" << endl;
        return -1;
    }

    string line;

    // Skip until POINTS
    while (getline(ob, line)) {
        if (line.find("POINTS") != string::npos)
            break;
    }

    // Read node coordinates
    for (int i = 0; i < nNODES; ++i) {
        ob >> X[i][0] >> X[i][1] >> X[i][2];
    }

    // Skip until CELLS
    while (getline(ob, line)) {
        if (line.find("CELLS") != string::npos)
            break;
    }

    // Read triangle connectivity
    for (int i = 0; i < nELEM; ++i) {
        int nodesPerCell;
        ob >> nodesPerCell;
        if (nodesPerCell != 3) {
            cerr << "Non-triangular cell found at line " << i + 1 << ": " << nodesPerCell << " nodes." << endl;
            return -1;
        }
        ob >> T[i][0] >> T[i][1] >> T[i][2];
    }

    // Skip until LOOKUP_TABLE to read nodal scalar data (U)
    while (getline(ob, line)) {
        if (line.find("LOOKUP_TABLE") != string::npos)
            break;
    }

    // Read nodal scalar solution
    for (int i = 0; i < nNODES; ++i) {
        ob >> U[i];
    }

    ob.close();

    // Sample output
    cout << "First node coordinates: " << X[0][0] << ", " << X[0][1] << ", " << X[0][2] << endl;
    cout << "First triangle: " << T[0][0] << ", " << T[0][1] << ", " << T[0][2] << endl;
    cout << "First nodal value U: " << U[0] << endl;

    //---------------------------------------
    //      Read velocities
    //---------------------------------------

    return 0;
}
