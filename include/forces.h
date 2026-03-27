#pragma once
#include "linalg.h"
#include "beam.h"
#include <string>
#include <vector>
#include <functional>

// A distributed load: q(x,t) in N/m, defined by a C-style function pointer
// Users express loads as compiled lambdas via the config
using LoadFn = std::function<double(double x, double t)>;

struct PointForce {
    std::function<double(double t)> x_fn;      // position along beam [m]
    std::function<double(double t)> mag_fn;     // magnitude [N]
    std::function<double(double t)> angle_fn;   // angle from axial [degrees]
};

struct ForceConfig {
    LoadFn qy;   // transverse distributed load [N/m]
    LoadFn qx;   // axial distributed load [N/m]
    std::vector<PointForce> point_forces;
};

// Assemble global force vector at time t
// Uses 4-point Gauss quadrature for distributed loads
Vec compute_force_vector(double t,
                         const ForceConfig& fc,
                         const AssembledSystem& sys,
                         const BeamConfig& cfg);

// Apply BCs: zero out constrained DOFs in K, M, F (penalty / elimination method)
void apply_bcs(Mat& K, Mat& M, Vec& F, const std::vector<int>& bc_dofs);
