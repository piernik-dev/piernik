# src/ – source tree map

The core Fortran modules live here; problem-specific modules (for example `problems/*/initproblem.F90`) live with their problems. Each file normally contains one module. The module name is *usually* but not always the file name
(`divB.F90`→`div_B`, `multigridmultipole.F90`→`multipole`, `signal_handler.F90`→`signalhandler`, `ppprofiling.F90`→`ppp`):
`setup` finds modules by scanning `module` statements, while message prefixes `[file_basename:routine]` follow the *file* name (`[divB:...]`). The directory structure is thematic only: `setup` flattens everything into
one `obj*/` directory, so there are no per-directory namespaces and module names must be globally unique.
Which files are built is decided by `! pulled by <SYMBOL>` comments (see [../python/AGENTS.md](../python/AGENTS.md)),
not by directory. Per-area details are in the `AGENTS.md` of each directory below.

## Layers (a module may `use` modules from its own or lower layers; avoid upward `use`)
1. `base/` (+ `base/mpi/`) – constants, units, MPI wrappers, global parameters, timers/profiling, containers.
2. `grid/` – domain, decomposition, grid containers (`cg`), lists/levels, AMR, load balancing, boundaries.
3. `fluids/`, `magnetic/` – fluid components (`fluidindex`, `flind`), their init, sources, cosmic rays, B-field tools.
4. `scheme/` – integration of the fluid equations (RTVD, HLLC, Riemann split/unsplit, CT).
5. `gravity/`, `multigrid/`, `particles/`, `shear/`, `fargo/`, `supernovae/`, `opencl/` – optional physics/solvers.
6. `IO/` – log/tsl output, HDF5 dumps and restarts (uses nearly everything, hence it sits on top).
Reality is not that clean (e.g. `base/timestep.F90`, `base/all_boundaries.F90` use fluids/grid; see "Proposed cleanups").

## Program flow (`base/piernik.F90`, `base/initpiernik.F90`)
`init_piernik` (MPI, namelists from `problem.par`, domain, grid, fluids, schemes, gravity, IO, problem's
`initproblem`) → loop: `time_step` → `fluid_update` (`scheme/fluidupdate.F90`) with sources, refinement
update, load balancing, `check_log/check_tsl/write_data` → `cleanup_piernik` (`base/finalizepiernik.F90`).
Users plug in through hooks declared in `base/user_hooks.F90` (set from `initproblem.F90`).

## Key concepts
- `cg` = grid container (a block of cells with guardcells), organised in lists (`cg_list_*`) and
  levels (`cg_level_*`); `leaves`, `all_cg` are the standard lists; loops are usually
  `cgl => leaves%first; do while (associated(cgl)) ... cg => cgl%cg ... cgl => cgl%nxt`.
- Fields are *named arrays* (`base/named_array*.F90`): `qna` (3D) and `wna` (4D) registries; the fluid state
  is `cg%u(:,:,:,:)` indexed via `flind` (`fluids/fluidindex.F90`), magnetic field `cg%b`.
- Preprocessor symbols (from `piernik.def`, derived in `base/piernik.h`): `IONIZED`, `NEUTRAL`, `DUST`,
  `MAGNETIC`, `COSM_RAYS`, `CRESP`, `GRAV`, `SELF_GRAV`, `MULTIGRID`, `NBODY`, `SHEAR`, `RESISTIVE`, `SN_SRC`,
  `THERM`, `PIERNIK_OPENCL`, `DEBUG`, ... full list: `doc/general/precompiler.F90`.
  Implications: `SELF_GRAV` ⇒ `GRAV` + `MULTIGRID`; `NBODY` ⇒ `SELF_GRAV`; `CRESP` ⇒ `COSM_RAYS`; `IONIZED` ⇒ `MAGNETIC`
  unless `NONMAGNETIC`.
- Run-time options are namelists (`doc/general/namelists.F90`), read in each module's `init_*`;
  a new namelist parameter needs: declaration, namelist group, MPI broadcast, validation.

## Small directories (no own AGENTS.md)
- `fargo/fargo.F90` – FARGO orbital advection for differentially rotating disks.
- `shear/shear.F90` – shearing-box boundaries (`SHEAR`).
- `supernovae/snsources.F90` – supernova insertion (`SN_SRC`).
- `opencl/piernikcl.F90` – OpenCL offload (`PIERNIK_OPENCL`, experimental).
- `utils/signal_handler.F90` – signal handling (Intel `ifport` based).
- `base/defines.c` is never compiled into the code: it only checks consistency of `piernik.def` (`#error` stops make).

## Proposed cleanups of the thematic split (not done; discuss before moving files)
Moving files is cheap for `setup` (it globs `src/**`) but changes `git blame`/history and docs, so batch moves.
1. `grid/` (51 files) mixes four topics: grid containers, lists/levels, **load balancing**
   (`cg_cost*`, `load_balance*`, `rebalance`, `ordering`, `sort_*`, `patch_list`) and **AMR**
   (`refinement.F90`, `refinement_flag`, `refinement_update` live in `grid/`, while URC criteria are in `grid/refinement/`).
   Suggest `grid/balance/` and moving the three `refinement*.F90` into `grid/refinement/`.
2. `base/` mixes generic infrastructure with physics/numerics: `timestep*`, `all_boundaries`, `fluidboundaries_funcs`,
   `coriolis`, `non_inertial`, `gridgeometry`, `mass_defect` depend on grid/fluids. Candidates for `grid/` or `scheme/`/`fluids/`.
   Also: profiling cluster (`ppp*`) → `base/ppp/`; containers (`cnt_array`, `log_bins`, `real_vector`,
   `sortable_list`, `primes`, `overlap`, `axes`) → `base/containers/`.
3. Multigrid is split in three places: generic in `multigrid/`, Poisson solver in `gravity/multigrid_*`, and
   diffusion solver in `fluids/cosmicrays/multigrid_diffusion.F90`. Suggest moving the `gravity/multigrid*`, `pcg`,
   `multipole_array`, `old_soln_list` files and the CR diffusion solver into `multigrid/`.
4. `scheme/ct.F90` (constrained transport) and `magnetic/` (divB, bfc_bcc) belong together; `magnetic/resist/` has one file.
5. `IO/hdf5/particles_io.F90` and `particles/particle_func.F90`, `particle_types.F90` are pulled by different
   symbols (`NBODY && HDF5`, `GRAV`) – verify the `pulled by` lines are intended.
6. Five single-file directories (`fargo`, `shear`, `supernovae`, `opencl`, `utils`) could be gathered (e.g. into
   `physics/` and `utils/`); single-file subdirs of `fluids/` are fine because they mirror fluid components.
7. Doc hygiene: modules lacking `\brief` (`piernik.F90`, `fluids_pub`, `cresp_variables`, `fluxlimiters`,
   `unsplit_*`), copy-pasted brief in `base/randomization.F90` (describes primes), inconsistent directory case
   (`IO`, `Riemann_Solver` vs lowercase elsewhere).
