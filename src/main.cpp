#include "beam.h"
#include "forces.h"
#include "solver.h"
#include "output.h"
#include <cmath>
#include <iostream>

// ============================================================
//  PROBLEM DEFINITION
//  Edit this section to set up your beam and loading.
//  All other modules (assembly, solver, output) stay untouched.
// ============================================================

static constexpr double PI = M_PI;

BeamConfig make_beam() {
    BeamConfig cfg;
    cfg.L      = 5.0;           // beam length [m]
    cfg.n_elem = 20;            // number of elements (more = more accurate but slower)
    cfg.E      = 200e9;         // Young's modulus [Pa]  (200 GPa = structural steel)
    cfg.I      = 8000e-8;       // 2nd moment of area [m^4]  (8000 cm^4)
    cfg.A      = 60e-4;         // cross-section area [m^2]  (60 cm^2)
    cfg.rho    = 7850.0;        // density [kg/m^3]
    cfg.zeta   = 0.02;          // damping ratio (2%)
    cfg.bc_left  = "fixed";     // "fixed" | "pinned" | "free"
    cfg.bc_right = "free";      // cantilever: fixed left, free right
    return cfg;
}

ForceConfig make_forces(const BeamConfig& cfg) {
    ForceConfig fc;

    // --- Distributed transverse load q_y(x, t) [N/m] ---
    // Oscillating uniform load: 1 kN/m * sin(2pi*t)
    fc.qy = [](double x, double t) -> double {
        (void)x;
        return 1000.0 * std::sin(2.0 * PI * t);
    };

    // --- Distributed axial load q_x(x, t) [N/m] ---
    // Zero by default
    fc.qx = [](double x, double t) -> double {
        (void)x; (void)t;
        return 0.0;
    };

    // --- Point forces ---
    // Example: transient point force at the free tip, directed at 80 deg from axial
    {
        PointForce pf;
        pf.x_fn     = [&cfg](double /*t*/) { return cfg.L; };         // at free end
        pf.mag_fn   = [](double t) { return 5000.0 * std::exp(-t); }; // decaying impulse
        pf.angle_fn = [](double /*t*/) { return 90.0; };              // pure transverse [deg]
        fc.point_forces.push_back(pf);
    }

    // Uncomment for a second moving force:
    // {
    //     PointForce pf;
    //     pf.x_fn     = [&cfg](double t) { return 0.5 * cfg.L * (1 + std::sin(t)); };
    //     pf.mag_fn   = [](double /*t*/) { return 3000.0; };
    //     pf.angle_fn = [](double t) { return 90.0 - 30.0 * std::sin(t); };
    //     fc.point_forces.push_back(pf);
    // }

    return fc;
}

SolverConfig make_solver_config() {
    SolverConfig sc;
    sc.t_start = 0.0;
    sc.t_end   = 2.0;   // [s]
    sc.n_steps = 500;   // time steps (more = smoother, slower)
    sc.beta    = 0.25;  // Newmark-beta (0.25 = unconditionally stable)
    sc.gamma   = 0.5;   // Newmark-gamma
    return sc;
}

// ============================================================
//  MAIN
// ============================================================
int main() {
    std::cout << "=== Beam FEA Transient Solver ===\n\n";

    BeamConfig   beam_cfg = make_beam();
    SolverConfig sol_cfg  = make_solver_config();

    std::cout << "[setup] L=" << beam_cfg.L << "m  n_elem=" << beam_cfg.n_elem
              << "  ndof=" << (beam_cfg.n_elem+1)*3 << "\n";
    std::cout << "[setup] BC: " << beam_cfg.bc_left << " | " << beam_cfg.bc_right << "\n";
    std::cout << "[setup] t=[" << sol_cfg.t_start << ", " << sol_cfg.t_end
              << "]  steps=" << sol_cfg.n_steps << "\n\n";

    AssembledSystem sys     = assemble(beam_cfg);
    ForceConfig     forces  = make_forces(beam_cfg);

    // Save every 2 steps -> ~250 frames for animation
    int save_every = std::max(1, sol_cfg.n_steps / 250);
    auto frames = newmark_solve(beam_cfg, sys, forces, sol_cfg, save_every);

    write_results(frames, beam_cfg, sys, "results");

    std::cout << "\n[done] Results written to results/\n";
    std::cout << "       Run:  python3 visualize.py\n";
    return 0;
}
