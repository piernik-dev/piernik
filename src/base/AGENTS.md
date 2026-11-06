# src/base – infrastructure

Lowest layer; used by everything. Parent: [../AGENTS.md](../AGENTS.md).

## Core definitions
- `constants.F90` – numeric constants, **enumerations** (solver ids, directions `xdim/ydim/zdim`, `PIERNIK_*` code-progress states, string lengths). Add new enum ids here.
- `units.F90` – physical constants and unit systems. `types.F90`, `axes.F90`, `problem_pub.F90`, `timestep_pub.F90`, and `../fluids/fluids_pub.F90` – small shared types/variables.
- `global.F90` – global run-time parameters (namelists `NUMERICAL_SETUP`, `PARALLEL_SETUP`, ...): `t`, `dt`, `nstep`, CFL, solver choice (`which_solver`: rtvd / hllc / riemann / unsplit), `divB_0_method`, limiters, etc.
- `piernik.h` – central header included by (almost) every file: pulls `piernik.def`, derives symbols, defines `RNG` (cell range macro), `MPIF`/`MPIFUN`. `defines.c` only validates `piernik.def` consistency.
- `piernik.F90` – main program; `initpiernik.F90` / `finalizepiernik.F90` – init and cleanup sequence (order matters: namelist-reading `init_*` calls depend on each other).

## Data structures
- `named_array.F90`, `named_array_list.F90` – 3D/4D named arrays and registries (`qna`, `wna`) used to attach fields to every grid container.
- Generic containers: `cnt_array`, `log_bins`, `real_vector`, `sortable_list` (extend it for sorting), `primes`, `overlap`, `randomization` (pulled by `RANDOMIZE`).
- `func.F90` – misc math helpers; `gridgeometry.F90` – Cartesian/cylindrical/polar geometry setup.

## Time stepping, boundaries, forces (depend on grid/fluids – see cleanup proposals in `../AGENTS.md`)
`timestep.F90` (collects dt limits from all modules), `repeatstep`, `all_boundaries.F90`, `fluidboundaries_funcs.F90`, `coriolis.F90` (`CORIOLIS`), `non_inertial.F90`, `mass_defect.F90`, `user_hooks.F90` (templates of user-overridable procedures: `problem_customize_solution`, `finalize_problem`, boundary hooks, ...).

## Diagnostics, profiling, resources
`diagnostics.F90`, `debug.F90` (`DEBUG`), `memory_usage.F90`, `timer.F90`, `procnames.F90`; PPP (precise parallel profiling): `ppprofiling.F90` (module `ppp`), `ppp_eventlist`, `ppp_events`, `mpi/ppp_mpi`.

## src/base/mpi – MPI layer
- `mpisetup.F90` – MPI init/cleanup and basic variables (`master`, rank, size); `mpi_wrappers.F90` – shared helpers.
- Typed wrappers with a uniform interface: `allreduce`, `bcast`, `barrier`, `isend_irecv`, `cmp_1D_mpi`.
- Request/tag bookkeeping for non-blocking exchanges: `req_array`, `tag_array(s)`, `sort_tags`, `wrapper_stats`.
- `MPIF` macro selects `mpi` or `mpi_f08` (`-DMPIF08`, auto-detected by the generated Makefile); keep code working with both.
  Files here contain `#warning` lines for Intel oneAPI; they are harmless (see `setup` notes).

## Gotchas
- Because `piernik.h` is included everywhere, touching it recompiles everything.
- Changing a module in `base` triggers a rebuild of most of the code; prefer small additions to `constants.F90` over new `use` chains.
- New modules here should not `use` `grid/`, `fluids/` etc.; the exceptions are the legacy ones listed above.
