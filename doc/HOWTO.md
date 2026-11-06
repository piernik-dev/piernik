# How to ... (recipes for frequent changes)

Each recipe lists *where to look* and *what must stay consistent*; copy an existing, similar case rather than writing from scratch
(`git log -S'<name>'` or `grep -rn` of a similar parameter shows a complete example). Verify the steps against the current code – it changes.
Source map: [../src/AGENTS.md](../src/AGENTS.md); vocabulary: [GLOSSARY.md](GLOSSARY.md).

After any recipe: `./setup <problem> -n`, a build, a short run of 1–2 standard problems (`sedov`, a problem using your feature), `make QA`.

## 1. Add a run-time parameter (namelist)
Pattern: see `init_cosmicrays` in `src/fluids/cosmicrays/initcosmicrays.F90` (or any `init_*` with `namelist /NAME/`).
1. Declare the variable in the module and add it to the `namelist /NAME/ ...` statement.
2. Set its **default** in the block before the `if (master)` part.
3. Add it to the transfer buffers in the `master` block (`rbuff`/`ibuff`/`lbuff`/`cbuff`) and the matching `slave` block that unpacks them. **Indices must match on both sides**, and the counters (`nl`, `nn`) must be updated (comments in the code say "must match the last ... index").
4. Document it next to the namelist (`!<` comment) – this is what Doxygen shows.
5. Only the master process reads `problem.par`; other ranks learn values through the buffers – a missing slave line gives *silently different values on different ranks*. Test with `mpirun -n 2`.
6. Parameter files in `problems/**/problem.par*` that set the old name will fail to parse if you rename/remove something (`grep -rn name problems`).

## 2. Add a new problem
See [../problems/AGENTS.md](../problems/AGENTS.md) ("Adding a problem"). Fast path: copy the closest problem, then edit `piernik.def` (physics), `initproblem.F90` (ICs), `problem.par`, `info`.

## 3. Add a user hook / custom behaviour for one problem
Hooks are procedure pointers in `src/base/user_hooks.F90` (e.g. `problem_customize_solution`, `problem_refine_derefine`, `late_initial_conditions`) and
`src/IO/dataio_user.F90` (`user_vars_hdf5`, `user_tsl`, `user_attrs_*`). Assign them in your `problem_pointers` (example: `problems/sedov/initproblem.F90`).
Contract: `doc/general/initproblem.F90`. A new hook = new pointer (initialised to `NULL()`) + an interface + a guarded call (`if (associated(...))`) at the right place in the main loop.

## 4. Add a field (array) to the grid
Fields are registered by name in the named-array registries: `call all_cg%reg_var(name, ...)` (examples: `src/gravity/gravity.F90` for `gpot`, with restart mode;
`src/gravity/pcg.F90`). Decide: does it need ghost-cell prolongation order (`ord_prolong`), is it *vital* (must be copied to new blocks at refinement), is it written to restart (`restart_mode`).
Access it through `cg%q(qna%ind(name))%arr` (3D) or `cg%w(wna%ind(name))%arr` (4D); cache the index at init time, not per cell.

## 5. Output a new variable to HDF5
Standard variable: add its name to the variable list handling in the `*_hdf5` module under `src/IO/` and keep names consistent with the Python tools in `python/` and `visual/`
(they look up datasets by name). Problem-specific variable: use the `user_vars_hdf5` hook (recipe 3) – no change in `src/` needed.
Time-series (`*_tsl`) columns: the `user_tsl` hook (`dataio_user::tsl_out`).

## 6. Add a boundary condition
Boundary types are enumerators `BND_*` in `src/base/constants.F90` and string names in the `bnd_xl ...` parameters of `BASE_DOMAIN` (`src/grid/domain.F90`).
Where they are applied: `src/grid/cg_list_bnd.F90`, `src/base/all_boundaries.F90`, `src/base/fluidboundaries_funcs.F90` (fluids) and the magnetic/gravity counterparts.
Search for an existing similar type (e.g. `BND_OUTH`) with `grep -rn BND_OUTH src` and add your type everywhere that one appears; unhandled types should `die`.

## 7. Add a source term / physics component
Sources are collected in `src/fluids/sources.F90`; physics components live in `src/fluids/<component>/` and are switched by a preprocessor symbol from `piernik.def`.
New files need a `! pulled by SYMBOL` comment so that `./setup` compiles them (see [../python/AGENTS.md](../python/AGENTS.md)); module names must be globally unique.
Check with `./setup <problem> -n` that the file is picked up only when intended. To verify that problems **without** your symbol still build, run `make allsetup`, then `make obj_A_* -j`.

## 8. Add a refinement criterion
User-defined: hook `problem_refine_derefine` or the URC framework (`src/grid/refinement/URC*.F90`; `URC_user`, `URC_Jeans`, `URC_nbody` are examples).
Framework criteria are configured through the `AMR` namelist (e.g. `refine_vars`). Sedov's problem-specific criterion uses `ref_thr` and `ref_eps` in `PROBLEM_CONTROL`; see `problems/sedov/initproblem.F90`.

## 9. Add a new compiler configuration
Copy `compilers/default.in`, adapt, test with `./setup sedov -c <name> -n` and a build. Do not edit `default.in` for a local setup
(see [../compilers/AGENTS.md](../compilers/AGENTS.md)).

## 10. Add a QA/CI check
Quick checks go to the `QA` group in the root `Makefile` (scripts in `bin/`, e.g. `checkmessages.sh`, `check_license_headers.sh`, `qa.py`).
Build/regression checks go to `artifacts`/`gold`; for gold tests see [../jenkins/AGENTS.md](../jenkins/AGENTS.md). Keep the target names stable.
