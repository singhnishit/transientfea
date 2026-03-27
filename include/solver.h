#pragma once
#include "linalg.h"
#include "beam.h"
#include "forces.h"
#include <vector>
#include <string>

struct SolverConfig {
    double t_start;
    double t_end;
    int    n_steps;
    double beta  = 0.25;   // Newmark-beta (0.25 = trapezoidal, unconditionally stable)
    double gamma = 0.5;    // Newmark-gamma
};

struct SimFrame {
    double t;
    Vec u;      // displacement [ndof]  (axial, transverse, rotation per node)
    Vec v;      // velocity     [ndof]
    Vec stress; // bending stress per element [n_elem]
};

// Stress recovery: axial strain + bending curvature at element midpoint
Vec compute_stress(const Vec& u, const AssembledSystem& sys,
                   const BeamConfig& cfg);

// Run full Newmark-beta time integration
// Returns one SimFrame per saved timestep (every `save_every` steps)
std::vector<SimFrame> newmark_solve(const BeamConfig& cfg,
                                    const AssembledSystem& sys,
                                    const ForceConfig& fc,
                                    const SolverConfig& sc,
                                    int save_every = 1);
