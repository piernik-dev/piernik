# src/scheme – integration of the fluid equations

Parent: [../AGENTS.md](../AGENTS.md).

## Driver (compiled always)
- `fluidupdate.F90` – `fluid_update`: advances the state by one step (pairs of sweeps in alternating directions for split schemes),
  calls sources, boundaries, CR/gravity/shear/FARGO hooks, and selects the solver by `which_solver`
  (`HLLC_SPLIT`, `RTVD_SPLIT`, `RIEMANN_SPLIT`, `RIEMANN_UNSPLIT` in `constants.F90`).
- `sweeps.F90` – one directional sweep over all blocks (with inter-block flux communication); `fluxtypes`; `timestep_retry`
  (repeat the step with smaller dt when CFL or positivity is violated).
- The solver is a **run-time** choice (`solver_str` in namelist `NUMERICAL_SETUP`), so all solver directories are built together.
  `RTVD` only supports CT for ∇·B = 0; Riemann solvers default to `HDC` (hyperbolic divergence cleaning).

## Solvers
| Directory | Scheme | Notes |
|-----------|--------|-------|
| `rtvd_split/` | relaxing TVD (`rtvd`, `solve_cg_rtvd`) | oldest, robust, 1D split; needs CT |
| `Riemann_Solver/` | MUSCL-type Riemann (split): `solve_cg_riemann`, `hlld` (Miyoshi–Kusano HLLD), `interpolations`, `fluxlimiters`, `hdc` | default |
| `unsplit_schemes/` | multi-D unsplit Riemann: `unsplit_fluidupdate`, `unsplit_sweeps`, `solve_cg_unsplit`, `unsplit_state_update`, `unsplit_source`, `unsplit_mag_modules` | newer; some files lack `\brief` |
| `hllc/` | HLLC split MUSCL-Hancock (`fluidupdate_hllc`) | **not maintained** (warning at run time) |
- `ct.F90` (`MAGNETIC`) – constrained transport of face-centred B (belongs conceptually with `../magnetic/`).

## Gotchas
- Block-level solvers (`solve_cg_*`) work on one `cg`; communication of fluxes at fine/coarse boundaries is in `grid/fc_fluxes`.
- After touching a solver run `make gold` and the `artifacts` tests, plus `sedov` and `otvortex` (MHD) by hand.
- Fluid variables are indexed through `flind` (`../fluids/fluidindex.F90`); never hard-code positions.
