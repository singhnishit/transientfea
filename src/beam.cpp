#include "beam.h"
#include <cmath>
#include <stdexcept>

BeamElement build_element(double E, double I, double A, double rho, double Le) {
    BeamElement el;
    el.Ke = mat_zeros(6, 6);
    el.Me = mat_zeros(6, 6);

    // Axial stiffness (DOFs 0 and 3)
    double EA = E * A / Le;
    el.Ke[0][0] =  EA; el.Ke[0][3] = -EA;
    el.Ke[3][0] = -EA; el.Ke[3][3] =  EA;

    // Bending stiffness (DOFs 1,2,4,5) 
    double EI = E * I;
    double k  = EI / (Le*Le*Le);
    // [v1, t1, v2, t2] -> global indices [1,2,4,5]
    int bi[4] = {1, 2, 4, 5};
    double kv[4][4] = {
        { 12,      6*Le,   -12,      6*Le},
        {  6*Le,   4*Le*Le, -6*Le,   2*Le*Le},
        {-12,     -6*Le,    12,     -6*Le},
        {  6*Le,   2*Le*Le, -6*Le,   4*Le*Le}
    };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            el.Ke[bi[i]][bi[j]] += k * kv[i][j];

    // Consistent mass matrix
    double rhoA = rho * A;

    // Axial mass (linear shape functions)
    el.Me[0][0] = rhoA * Le / 3.0;
    el.Me[0][3] = rhoA * Le / 6.0;
    el.Me[3][0] = rhoA * Le / 6.0;
    el.Me[3][3] = rhoA * Le / 3.0;

    // Bending mass (Hermite cubic shape functions)
    double m = rhoA * Le / 420.0;
    double mv[4][4] = {
        { 156,      22*Le,     54,     -13*Le},
        {  22*Le,    4*Le*Le,  13*Le,   -3*Le*Le},
        {  54,      13*Le,    156,     -22*Le},
        { -13*Le,   -3*Le*Le, -22*Le,    4*Le*Le}
    };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            el.Me[bi[i]][bi[j]] += m * mv[i][j];

    return el;
}

static void apply_node_bc(std::vector<int>& bc_dofs, int node, const std::string& type) {
    int base = node * 3;
    if (type == "fixed") {
        bc_dofs.push_back(base);
        bc_dofs.push_back(base + 1);
        bc_dofs.push_back(base + 2);
    } else if (type == "pinned") {
        bc_dofs.push_back(base);
        bc_dofs.push_back(base + 1);
    }
    // "free" -> no constraints
}

AssembledSystem assemble(const BeamConfig& cfg) {
    int n  = cfg.n_elem;
    int nn = n + 1;
    int ndof = nn * 3;
    double Le = cfg.L / n;

    Mat K = mat_zeros(ndof, ndof);
    Mat M = mat_zeros(ndof, ndof);

    for (int e = 0; e < n; e++) {
        BeamElement el = build_element(cfg.E, cfg.I, cfg.A, cfg.rho, Le);
        int dofs[6] = {e*3, e*3+1, e*3+2, (e+1)*3, (e+1)*3+1, (e+1)*3+2};
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) {
                K[dofs[i]][dofs[j]] += el.Ke[i][j];
                M[dofs[i]][dofs[j]] += el.Me[i][j];
            }
    }

    std::vector<int> bc_dofs;
    apply_node_bc(bc_dofs, 0,    cfg.bc_left);
    apply_node_bc(bc_dofs, nn-1, cfg.bc_right);

    AssembledSystem sys;
    sys.K       = K;
    sys.M       = M;
    sys.ndof    = ndof;
    sys.n_nodes = nn;
    sys.Le      = Le;
    sys.bc_dofs = bc_dofs;
    return sys;
}
