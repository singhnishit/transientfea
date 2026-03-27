#pragma once
#include "linalg.h"

// Euler-Bernoulli beam element
// DOFs per node: [u (axial), v (transverse), theta (rotation)]
// Element DOF order: [u1, v1, t1, u2, v2, t2]

struct BeamElement {
    Mat Ke;  // 6x6 stiffness
    Mat Me;  // 6x6 consistent mass
};

// Build a single element given material/section properties and element length
BeamElement build_element(double E, double I, double A, double rho, double Le);

struct BeamConfig {
    double L;       // total length [m]
    int n_elem;     // number of elements
    double E;       // Young's modulus [Pa]
    double I;       // second moment of area [m^4]
    double A;       // cross-section area [m^2]
    double rho;     // density [kg/m^3]
    double zeta;    // damping ratio

    // boundary conditions: "fixed", "pinned", "free"
    std::string bc_left;
    std::string bc_right;
};

struct AssembledSystem {
    Mat K;   // global stiffness [ndof x ndof]
    Mat M;   // global mass     [ndof x ndof]
    int ndof;
    int n_nodes;
    double Le;  // element length
    std::vector<int> bc_dofs;  // constrained DOF indices
};

AssembledSystem assemble(const BeamConfig& cfg);
