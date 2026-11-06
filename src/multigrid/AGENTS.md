# src/multigrid – multigrid infrastructure

Parent: [../AGENTS.md](../AGENTS.md). Pulled by `MULTIGRID` (set automatically by `SELF_GRAV`); FFT use needs FFTW (`NO_FFT` disables it).
Users of this infrastructure: Poisson solver in `../gravity/multigrid_*`, implicit CR diffusion in `../fluids/cosmicrays/multigrid_diffusion.F90`.

- `multigrid.F90` – init/cleanup, namelist `MULTIGRID_SOLVER` (`level_depth`, prolongation order, V-cycle verbosity), allocates levels/fields.
- `multigridvars.F90`, `multigrid_helpers.F90`, `multigrid_vstats.F90` – shared variables/types, helpers, V-cycle convergence statistics.
- Laplace operators (selected in `multigrid_Laplace.F90`): `multigrid_Laplace2` (2nd order), `multigrid_Laplace4` (4th, 5-point wide),
  `multigrid_Laplace4M` (4th order compact, Mehrstellen).
- Grid data for the coarse levels live in the grid containers (see `grid/grid_container_ext.F90`).

Tests: `maclaurin` (accuracy of the Poisson solver), `jeans`. Convergence regressions are easy to miss – compare residual histories in the `.log` file.
