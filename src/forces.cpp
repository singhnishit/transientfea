#include "forces.h"
#include <cmath>
#include <algorithm>

Vec compute_force_vector(double t,
                         const ForceConfig& fc,
                         const AssembledSystem& sys,
                         const BeamConfig& cfg) {
    int ndof = sys.ndof;
    int n    = cfg.n_elem;
    double Le = sys.Le;
    Vec F = vec_zeros(ndof);

    // --- Distributed loads via 4-point Gauss quadrature ---
    // Gauss points on [-1,1], mapped to [0,1] parametric space
    const int NGP = 4;
    const double gp[NGP] = {0.0694318, 0.3300095, 0.6699905, 0.9305682};
    const double gw[NGP] = {0.1739274, 0.3260726, 0.3260726, 0.1739274};

    for (int e = 0; e < n; e++) {
        double x1 = e * Le;
        int dofs[6] = {e*3, e*3+1, e*3+2, (e+1)*3, (e+1)*3+1, (e+1)*3+2};

        for (int g = 0; g < NGP; g++) {
            double xi = gp[g];            // parametric coord in [0,1]
            double xg = x1 + xi * Le;    // physical coord
            double w  = gw[g] * Le;      // weight * jacobian

            double qy_val = fc.qy ? fc.qy(xg, t) : 0.0;
            double qx_val = fc.qx ? fc.qx(xg, t) : 0.0;

            // Linear shape functions for axial (u)
            double N1_ax = 1.0 - xi;
            double N2_ax = xi;

            // Hermite shape functions for transverse (v)
            // N1 = 1 - 3xi^2 + 2xi^3
            // N2 = xi*Le*(1 - 2xi + xi^2)  -- actually: Le*(xi - 2xi^2 + xi^3)
            // N3 = 3xi^2 - 2xi^3
            // N4 = xi*Le*(xi^2 - xi)        -- Le*(xi^2 - xi^3) ... sign convention below
            double xi2 = xi*xi, xi3 = xi2*xi;
            double H1 =  1.0 - 3.0*xi2 + 2.0*xi3;
            double H2 =  Le  * (xi - 2.0*xi2 + xi3);
            double H3 =  3.0*xi2 - 2.0*xi3;
            double H4 =  Le  * (xi3 - xi2);

            // Axial DOFs (0 and 3)
            F[dofs[0]] += qx_val * N1_ax * w;
            F[dofs[3]] += qx_val * N2_ax * w;

            // Transverse DOFs (1,2,4,5)
            F[dofs[1]] += qy_val * H1 * w;
            F[dofs[2]] += qy_val * H2 * w;
            F[dofs[4]] += qy_val * H3 * w;
            F[dofs[5]] += qy_val * H4 * w;
        }
    }

    // --- Point forces: decompose into axial and transverse ---
    int n_nodes = sys.n_nodes;
    for (const auto& pf : fc.point_forces) {
        double xpos  = pf.x_fn(t);
        double mag   = pf.mag_fn(t);
        double angle = pf.angle_fn(t) * M_PI / 180.0; // deg -> rad

        double Fx = mag * std::cos(angle);  // axial component
        double Fy = mag * std::sin(angle);  // transverse component

        // Find nearest node
        int node = (int)std::round(xpos / Le);
        node = std::max(0, std::min(n_nodes - 1, node));

        F[node*3]     += Fx;
        F[node*3 + 1] += Fy;
    }

    return F;
}

void apply_bcs(Mat& K, Mat& M, Vec& F, const std::vector<int>& bc_dofs) {
    int ndof = K.size();
    for (int d : bc_dofs) {
        for (int i = 0; i < ndof; i++) {
            K[d][i] = 0.0; K[i][d] = 0.0;
            M[d][i] = 0.0; M[i][d] = 0.0;
        }
        K[d][d] = 1.0;
        M[d][d] = 1.0;
        F[d]    = 0.0;
    }
}
