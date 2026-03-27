#include "output.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

static std::ofstream open_file(const std::string& path) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open output file: " + path);
    f << std::scientific << std::setprecision(8);
    return f;
}

void write_results(const std::vector<SimFrame>& frames,
                   const BeamConfig& cfg,
                   const AssembledSystem& sys,
                   const std::string& out_dir) {
    int n_nodes = sys.n_nodes;
    int n_elem  = cfg.n_elem;
    int n_frames = frames.size();

    // ---- deflection.csv ----
    // rows = frames, columns = nodes: t, v0, v1, ..., v_{n_nodes-1}
    {
        auto f = open_file(out_dir + "/deflection.csv");
        f << "t";
        for (int i = 0; i < n_nodes; i++) f << ",node_" << i;
        f << "\n";
        for (const auto& fr : frames) {
            f << fr.t;
            for (int i = 0; i < n_nodes; i++) f << "," << fr.u[i*3 + 1];
            f << "\n";
        }
        std::cout << "[output] wrote " << out_dir << "/deflection.csv\n";
    }

    // ---- axial.csv ----
    {
        auto f = open_file(out_dir + "/axial.csv");
        f << "t";
        for (int i = 0; i < n_nodes; i++) f << ",node_" << i;
        f << "\n";
        for (const auto& fr : frames) {
            f << fr.t;
            for (int i = 0; i < n_nodes; i++) f << "," << fr.u[i*3];
            f << "\n";
        }
        std::cout << "[output] wrote " << out_dir << "/axial.csv\n";
    }

    // ---- stress.csv ----
    {
        auto f = open_file(out_dir + "/stress.csv");
        f << "t";
        for (int i = 0; i < n_elem; i++) f << ",elem_" << i;
        f << "\n";
        for (const auto& fr : frames) {
            f << fr.t;
            for (int i = 0; i < n_elem; i++) f << "," << fr.stress[i];
            f << "\n";
        }
        std::cout << "[output] wrote " << out_dir << "/stress.csv\n";
    }

    // ---- metadata.csv ----
    {
        auto f = open_file(out_dir + "/metadata.csv");
        f << std::defaultfloat;
        f << "key,value\n";
        f << "L,"       << cfg.L       << "\n";
        f << "n_elem,"  << cfg.n_elem  << "\n";
        f << "n_nodes," << n_nodes     << "\n";
        f << "n_frames,"<< n_frames    << "\n";
        f << "Le,"      << sys.Le      << "\n";
        f << "E,"       << cfg.E       << "\n";
        f << "I,"       << cfg.I       << "\n";
        f << "A,"       << cfg.A       << "\n";
        f << "rho,"     << cfg.rho     << "\n";
        f << "t_start," << frames.front().t << "\n";
        f << "t_end,"   << frames.back().t  << "\n";
        std::cout << "[output] wrote " << out_dir << "/metadata.csv\n";
    }
}
