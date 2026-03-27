# Finite Element Analysis with transient forces

An implementation of the finite element method for 3D meshes with transient loads. The project is undertaken under [Dr. Gaurav Sharma](https://www.bits-pilani.ac.in/pilani/gaurav-sharma/) and will be integrated with a fluid solver for a fluid-structure interaction engine upon completion.

## Compiling and running
To run the code : 
```bash
unzip beam_fea_solver.zip && cd beam_fea
make visualize          # compiles + runs solver + renders animation
# or step by step:
make                    # compile only
make run                # compile + run solver -> results/*.csv
python3 visualize.py    # render animation -> results/beam_animation.gif
```
## Customising the problem
### To edit the beam properties
Navigate to `make_beam()` in `src/main.cpp` :
```cpp
BeamConfig make_beam() {
    BeamConfig cfg;
    cfg.L      = 5.0;       // total beam length in metres
    cfg.n_elem = 20;        // number of finite elements — more = more accurate, slower
    cfg.E      = 200e9;     // Young's modulus in Pascals  (200e9 = 200 GPa = structural steel)
    cfg.I      = 8000e-8;   // second moment of area in m^4  (8000e-8 = 8000 cm^4)
    cfg.A      = 60e-4;     // cross-section area in m^2  (60e-4 = 60 cm^2)
    cfg.rho    = 7850.0;    // density in kg/m^3
    cfg.zeta   = 0.02;      // damping ratio  (0.02 = 2%, typical for steel)
    cfg.bc_left  = "fixed"; // left end boundary condition
    cfg.bc_right = "free";  // right end boundary condition
}
```
Where boundary conditions are `fixed`, `free` and `pinned`. 
### To edit point forces
Navigate to the `PointForce` entries in `make_forces()` in `src/main.cpp` :
```cpp
{
    PointForce pf;
    pf.x_fn     = [&cfg](double t) { return cfg.L; };          // position along beam [m]
    pf.mag_fn   = [](double t)     { return 5000.0; };         // magnitude [N]
    pf.angle_fn = [](double t)     { return 90.0; };           // angle from axial axis [degrees]
    fc.point_forces.push_back(pf);
}
```
### To edit distributed loads
Navigate to `fc.qy` (transverse loading) and `fc.qx` (axial loading) in `make_forces()` :
```cpp
fc.qy = [](double x, double t) -> double {
    return 1000.0 * std::sin(2.0 * PI * t);  // your expression here
};

fc.qx = [](double x, double t) -> double {
    return 0.0;  // your expression here
};
```
### To edit time parameters
Navigate to `make_solver_config()` :
```cpp
SolverConfig make_solver_config() {
    SolverConfig sc;
    sc.t_start = 0.0;   // start time [s]
    sc.t_end   = 2.0;   // end time [s]
    sc.n_steps = 500;   // number of time steps
}
```

## Example config : 
### Point force(s)
```cpp
// --- position ---
pf.x_fn = [&cfg](double t) { return cfg.L; };               // fixed at free tip
pf.x_fn = [&cfg](double t) { return cfg.L / 2.0; };         // fixed at midspan
pf.x_fn = [&cfg](double t) { return 0.5 * cfg.L * (1.0 + std::sin(t)); }; // moves along beam

// --- magnitude ---
pf.mag_fn = [](double t) { return 5000.0; };                          // constant
pf.mag_fn = [](double t) { return 5000.0 * std::sin(2.0 * PI * t); }; // oscillating
pf.mag_fn = [](double t) { return 10000.0 * std::exp(-2.0 * t); };    // decaying impulse
pf.mag_fn = [](double t) { return t < 0.5 ? 8000.0 : 0.0; };         // step load, on for 0.5s then off

// --- angle ---
pf.angle_fn = [](double t) { return 90.0; };                          // always transverse
pf.angle_fn = [](double t) { return 0.0; };                           // always axial
pf.angle_fn = [](double t) { return 45.0; };                          // 45 degrees
pf.angle_fn = [](double t) { return 90.0 - 45.0 * std::sin(t); };    // rotating over time
```
### Distributed force(s)
```cpp
fc.qy = [](double x, double t) -> double {
    return 1000.0 * std::sin(2.0 * PI * t);  // your expression here
};

fc.qx = [](double x, double t) -> double {
    return 0.0;  // your expression here
};
```

