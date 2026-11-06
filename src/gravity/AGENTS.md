# src/gravity – gravitational potential, self-gravity solvers, hydrostatic equilibrium

Parent: [../AGENTS.md](../AGENTS.md). Pulled by `GRAV` (external potential) or `MULTIGRID && SELF_GRAV` (Poisson solver).
Generic multigrid machinery is in [../multigrid/](../multigrid/AGENTS.md).

- `gravity.F90` – namelist `GRAVITY`; external potentials (`external_gp`: null, linear, uniform, point mass variants, Roche, user),
  gravitational acceleration, source term for the fluids, `grav_pot_3d` computes the 3D potential.
- `hydrostatic.F90` – builds vertical hydrostatic equilibrium (`hydrostatic_zeq*`, `HYDROSTATIC_V2`).
- **Poisson solver** (`SELF_GRAV` ⇒ `MULTIGRID`): `multigrid_gravity` (V-cycle driver, namelist `MULTIGRID_GRAVITY`, boundary conditions),
  `multigrid_gravity_helper`, `multigrid_old_soln` + `old_soln_list` (history of solutions for initial guesses; affects restart format),
  `pcg` (preconditioned CG on the coarsest level).
- Isolated boundaries: `multigrid_monopole` (module `mg_monopole`), `multigridmultipole` (module `multipole`) with `multipole_array`.
- Particles couple through `particles/particle_gravity` (needs `NBODY` ⇒ `SELF_GRAV`).

Tests: `maclaurin`, `jeans`, `selfgrav_clump`, `lingrav` problems.
Gotcha: file name ≠ module name for the multipole files – grep `^ *module` before searching for `use` lines.
