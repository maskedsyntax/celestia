# Celestia: N-Body Gravitational Simulation

A Fortran-based N-Body simulation implementing Newtonian gravity with Barnes-Hut optimization for O(N log N) performance.

## Core Features
- **Barnes-Hut Algorithm**: Utilizes a 3D Octree to approximate long-range gravitational forces.
- **Velocity Verlet Integration**: Second-order symplectic integrator for stable energy conservation.
- **OpenMP Parallelism**: Parallelized force calculation and tree traversal.
- **Collision Handling**: Inelastic merging of bodies based on conservation of mass and momentum.
- **Softening Length**: Implemented to prevent numerical singularities during close encounters.

## Requirements
- `gfortran` (with OpenMP support)
- `make`
- `python3` (for data conversion only)

## Build and Run
1. Compile the simulation:
   ```bash
   make
   ```
2. Execute the simulation:
   ```bash
   make run
   ```
   The binary generates `galaxy_step.csv` containing state data (time, id, mass, pos, vel).

## Visualization
The project includes a zero-dependency web viewer.
1. Convert the simulation output to JavaScript format:
   ```bash
   python3 to_js.py
   ```
2. Open `viewer.html` in any web browser.
   - **Space**: Play/Pause
   - **Auto-scaling**: The viewer dynamically adjusts zoom based on particle distribution.
