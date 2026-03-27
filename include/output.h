#pragma once
#include "solver.h"
#include "beam.h"
#include <string>
#include <vector>

// Write all frames to CSV files in the given directory.
// Produces:
//   results/deflection.csv   -- transverse displacement v per node per frame
//   results/axial.csv        -- axial displacement u per node per frame
//   results/stress.csv       -- bending stress per element per frame
//   results/metadata.csv     -- beam geometry and time info for the plotter
void write_results(const std::vector<SimFrame>& frames,
                   const BeamConfig& cfg,
                   const AssembledSystem& sys,
                   const std::string& out_dir);
