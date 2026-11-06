# src/particles – N-body particles (`NBODY`)

Parent: [../AGENTS.md](../AGENTS.md). `NBODY` ⇒ `SELF_GRAV` ⇒ `MULTIGRID`; particles are gravitating tracers (stars, dark matter) living in grid containers.
Namelist `PARTICLES` is read in `particle_pub`.

- `particle_pub.F90` – main module: init, namelist (integrator, interpolation scheme, ...), public state.
- `particle_types.F90`, `particle_utils.F90` – particle/list types, creation, migration between blocks/ranks.
- `particle_solvers.F90` – time integrators (`leapfrog_2ord`, `hermit_4ord`); `particle_timestep.F90` – dt limits; `particle_func.F90` – helpers.
- `particle_gravity.F90` – mass deposit onto the grid (density for the Poisson solver) and acceleration back (via `particle_maps`).
- `particle_maps.F90` – grid ↔ particle interpolation (see `interpolation_scheme`).
- `star_formation.F90` – star formation and feedback; `particle_diag.F90` – diagnostics.
- Related: `grid/refinement/URC_nbody.F90` (refine by particle count), `IO/hdf5/particles_io.F90`.

Tests: `3body_orig`, `particle_map_test`. Note some files are `pulled by GRAV` rather than `NBODY` (`particle_types`, `particle_func`) – see cleanup list in `../AGENTS.md`.
