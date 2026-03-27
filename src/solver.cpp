#include "solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

Vec compute_stress(const Vec& u, const AssembledSystem& sys,
                   const BeamConfig& cfg) {
    int n  = cfg.n_elem;
    double Le = sys.Le;
    // radius of gyration used as extreme fibre distance estimate
    double y_max = std::sqrt(cfg.I / cfg.A);

    Vec stress(n);
    for (int e = 0; e < n; e++) {
        double u1 = u[e*3],     v1 = u[e*3+1], t1 = u[e*3+2];
        double u2 = u[(e+1)*3], v2 = u[(e+1)*3+1], t2 = u[(e+1)*3+2];

        // Axial strain at element midpoint (constant along element)
        double eps_axial = (u2 - u1) / Le;

        // Bending curvature at xi = 0.5 (midpoint)
        // kappa = d^2v/dx^2 from Hermite shape functions
        double xi = 0.5;
        // Second derivatives of Hermite functions w.r.t. x, evaluated at xi
        // d2H1/dx2 = (1/Le^2)*(- 6 + 12*xi)
        // d2H2/dx2 = (1/Le)  *(- 4 + 6*xi)
        // d2H3/dx2 = (1/Le^2)*(  6 - 12*xi)
        // d2H4/dx2 = (1/Le)  *(- 2 + 6*xi)
        double Le2 = Le * Le;
        double d2H1 = (-6.0 + 12.0*xi) / Le2;
        double d2H2 = (-4.0 +  6.0*xi) / Le;
        double d2H3 = ( 6.0 - 12.0*xi) / Le2;
        double d2H4 = (-2.0 +  6.0*xi) / Le;

        double kappa = d2H1*v1 + d2H2*t1 + d2H3*v2 + d2H4*t2;

        // Stress = E*(eps_axial + y_max * kappa)
        stress[e] = cfg.E * (eps_axial + y_max * kappa);
    }
    return stress;
}

std::vector<SimFrame> newmark_solve(const BeamConfig& cfg,
                                    const AssembledSystem& sys,
                                    const ForceConfig& fc,
                                    const SolverConfig& sc,
                                    int save_every) {
    int ndof  = sys.ndof;
    double dt = (sc.t_end - sc.t_start) / sc.n_steps;
    double b  = sc.beta, g = sc.gamma;

    std::cout << "[solver] Newmark-beta: beta=" << b << " gamma=" << g
              << "  dt=" << dt << "  steps=" << sc.n_steps << "\n";

    // Rayleigh damping: C = alpha_M * M  (mass-proportional)
    // Estimate first natural freq from diagonal to get alpha_M
    // omega1 ~ sqrt(K[1][1] / M[1][1]) -- rough estimate
    double K11 = sys.K[1][1], M11 = sys.M[1][1];
    double omega1 = (M11 > 1e-20) ? std::sqrt(K11 / M11) : 1.0;
    double alpha_M = 2.0 * cfg.zeta * omega1;
    std::cout << "[solver] Rayleigh mass damping alpha_M=" << alpha_M
              << " (omega1~" << omega1 << " rad/s)\n";

    // Effective stiffness: K_eff = K + (1/(b*dt^2)) * M + (g/(b*dt)) * C
    // With C = alpha_M * M:
    // K_eff = K + [1/(b*dt^2) + g*alpha_M/(b*dt)] * M
    double coeff_M = 1.0/(b*dt*dt) + g*alpha_M/(b*dt);
    Mat Keff = mat_zeros(ndof, ndof);
    for (int i = 0; i < ndof; i++)
        for (int j = 0; j < ndof; j++)
            Keff[i][j] = sys.K[i][j] + coeff_M * sys.M[i][j];

    // Apply BCs to Keff (only structural DOFs change; BC rows/cols set to identity)
    for (int d : sys.bc_dofs) {
        for (int i = 0; i < ndof; i++) { Keff[d][i] = 0.0; Keff[i][d] = 0.0; }
        Keff[d][d] = 1.0;
    }

    // Initial conditions
    Vec u = vec_zeros(ndof);
    Vec v = vec_zeros(ndof);
    Vec a = vec_zeros(ndof);

    // Compute initial acceleration from F(t=0)
    {
        Vec F0 = compute_force_vector(sc.t_start, fc, sys, cfg);
        // Apply BCs to F0
        for (int d : sys.bc_dofs) F0[d] = 0.0;
        // M * a0 = F0 - K*u0 - C*v0 = F0 (since u0=v0=0)
        // Approximate: a0[i] = F0[i] / M[i][i]  (diagonal mass approx for init)
        for (int i = 0; i < ndof; i++) {
            double mii = sys.M[i][i];
            a[i] = (mii > 1e-20) ? F0[i] / mii : 0.0;
        }
        for (int d : sys.bc_dofs) a[d] = 0.0;
    }

    std::vector<SimFrame> frames;
    int total_frames = sc.n_steps / save_every + 1;
    frames.reserve(total_frames);

    // Save frame 0
    SimFrame f0;
    f0.t = sc.t_start; f0.u = u; f0.v = v;
    f0.stress = compute_stress(u, sys, cfg);
    frames.push_back(f0);

    int progress_interval = std::max(1, sc.n_steps / 20);

    for (int step = 1; step <= sc.n_steps; step++) {
        double t_new = sc.t_start + step * dt;

        // Predictors
        Vec u_pred = vec_add(vec_add(u, vec_scale(v, dt)),
                             vec_scale(a, dt*dt*(0.5 - b)));
        Vec v_pred = vec_add(v, vec_scale(a, dt*(1.0 - g)));

        // External force at t_new
        Vec F_ext = compute_force_vector(t_new, fc, sys, cfg);
        for (int d : sys.bc_dofs) F_ext[d] = 0.0;

        // Effective RHS: F_eff = F_ext + M*(1/(b*dt^2))*u_pred + ...
        // Full form: F_eff = F_ext - K*u_pred - C*v_pred + M*(1/(b*dt^2))*delta_u_pred
        // Equivalently, solve for delta_u = u_new - u_pred:
        // K_eff * delta_u = F_ext - K*u_pred - C*v_pred
        // Then u_new = u_pred + delta_u

        Vec Ku_pred = mat_mul_vec(sys.K, u_pred);
        // C*v_pred = alpha_M * M * v_pred
        Vec Cv_pred = vec_scale(mat_mul_vec(sys.M, v_pred), alpha_M);
        Vec RHS(ndof);
        for (int i = 0; i < ndof; i++)
            RHS[i] = F_ext[i] - Ku_pred[i] - Cv_pred[i];
        for (int d : sys.bc_dofs) RHS[d] = 0.0;

        Vec du = lu_solve(Keff, RHS);
        for (int d : sys.bc_dofs) du[d] = 0.0;

        Vec u_new = vec_add(u_pred, du);
        Vec a_new = vec_scale(du, 1.0/(b*dt*dt));
        Vec v_new(ndof);
        for (int i = 0; i < ndof; i++)
            v_new[i] = v_pred[i] + g*dt*a_new[i];

        u = u_new; v = v_new; a = a_new;

        if (step % save_every == 0) {
            SimFrame fr;
            fr.t = t_new; fr.u = u; fr.v = v;
            fr.stress = compute_stress(u, sys, cfg);
            frames.push_back(fr);
        }

        if (step % progress_interval == 0) {
            std::cout << "[solver] step " << step << "/" << sc.n_steps
                      << "  t=" << t_new << "\n";
        }
    }

    std::cout << "[solver] Done. " << frames.size() << " frames saved.\n";
    return frames;
}
