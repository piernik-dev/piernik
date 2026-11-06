# problems/ – computational problems (initial conditions + configuration)

Parent: [../AGENTS.md](../AGENTS.md). ~40 top-level problems, ~190 parameter files. A problem is selected by `./setup <path-relative-to-problems>`
(e.g. `sedov`, `sedov/MHD_blast_wave`, `2body/3body`, `testing_and_debuging/chimaera`). Mechanics of collecting the files: [../python/AGENTS.md](../python/AGENTS.md).

## Anatomy of a problem directory
| File | Role |
|------|------|
| `piernik.def` | **compile-time** choices: `#define IONIZED`, `GRAV`, `SELF_GRAV`, `COSM_RAYS`, ... (symbols documented in `doc/general/precompiler.F90`; derived ones in `src/base/piernik.h`) |
| `initproblem.F90` | module `initproblem` with `read_problem_par` (namelist `PROBLEM_CONTROL`), `problem_initial_conditions`, `problem_pointers` (sets user hooks from `src/base/user_hooks.F90` / `src/IO/dataio_user.F90`); contract in `doc/general/initproblem.F90` |
| `problem.par` | **run-time** parameters, one namelist per section (`BASE_DOMAIN`, `MPI_BLOCKS`, `AMR`, `OUTPUT_CONTROL`, `NUMERICAL_SETUP`, ...). Variants: `problem.par.<variant>` selected by `./setup -p problem.par.<variant>`. Common ones: `.build` (cheap, for compile/CI tests), `.gold` (gold tests), `.AMR`, `.restart_test_v2` |
| `info` | one-line description shown by `./setup --problems` (missing for some) |
| `Makefile.in` | optional extra make rules; **not inherited** by subproblems |
| `*.py`, `*.sh`, `*.gnuplot`, `README.md`, `artifacts.md` | problem-specific analysis helpers/notes (e.g. `particle_error.py` used by `make 3body`, `jeans.gnuplot` by `make Jeans`) |

**Inheritance:** a subdirectory (`sedov/MHD_blast_wave`) takes its own files first, then fills in missing ones from the parent directory (except `info`, `Makefile.in`).
`overrides` (a comment/directive in a problem's `.F90`) lets a problem replace a file from `src/` of the same basename; see `override_test/`. Use it only for experiments.

## Problems used by tests
| Test | Problem |
|------|---------|
| artifacts: `Jeans`, `Maclaurin`, `CPAW`, `IOv2` (`advection_test`), `3body` (`2body/3body`), `noHDF5` (`crtest`), `dep`/`I64` (`testing_and_debuging/chimaera`) | `jeans`, `maclaurin`, `cpaw`, `advection_test`, ... |
| gold (`jenkins/gold_configs`) | `mcrtest` (`mcrtest_new`), `mcrtest/mcrtest_cresp`, `mcrwind`, `sedov/MHD_blast_wave` (AMR), `tearing` (`resist`), `streaming_global` |
Parameter files used by tests must keep working: changing `problem.par.gold`/`.build` or initial conditions changes gold results (see [../jenkins/AGENTS.md](../jenkins/AGENTS.md)).

## Selected problems (see each `info`)
`sedov` (Sedov–Taylor; also standard performance benchmark), `maclaurin`/`jeans`/`selfgrav_clump`/`lingrav` (self-gravity), `otvortex`/`cpaw`/`field_loop_mhd`/`magnetic_rotor`/`khi`/`blob`/`tearing`/`shock1d_rj`/`sod2d`/`2D_RIEMANN` (hydro/MHD classics),
`crtest`/`mcrtest`/`mcrwind` (cosmic rays, CRESP), `dust*`/`stream*`/`streaming_global` (dust, streaming instability), `kepler`/`roche`/`galdiskmao4`/`wt4`/`fargo_test` (discs, rotation), `2body`/`particle_map_test` (N-body),
`therm_inst`, `turbulence`, `wind`, `slab`, `divvtest`, `advection_test`.
`testing_and_debuging/` is a drawer of old developer experiments ("not guaranteed to be useful or work anymore"); its `chimaera` includes nearly all of `src/` and is the *compile-everything* test. It has `.skipauto` so `make allsetup` skips it.

## Adding a problem
Copy the closest existing one; provide `piernik.def`, `initproblem.F90` (start from `doc/general/new_source_file.F90` for the licence header), `problem.par`, `info`.
Keep the default `problem.par` cheap (it is what users try first); add `problem.par.build` if a CI compile check is meaningful.
Verify: `./setup <name> -n`, a build and a short run, `make QA` (the checks also lint `problems/**/*.F90`). `make allsetup` sets up every problem without `.skipauto`.

## Observations
- Many problems have old Python helpers using `yt.mods`/`from pylab import *`/matplotlib `cairo` backend (8 files, e.g. `mcrtest/piernik_problem.py`) – they fail with current yt/matplotlib. Port or delete.
- `piernik_problem.py` exists in 4 problems but is not hooked to any test.
- Problem names are inconsistent (`testing_and_debuging` is misspelled; `stream` vs `stream3D` vs `streaming_global`; `2body` contains `3body` and `testing_and_debuging/N-body` exists as well). Renaming breaks `.setup.call` files, gold configs and user scripts – do it together.
- Some problems have `powerbug/burn.sh` / `fit_res_cfl.sh` scripts with hard-coded assumptions – read before running.
- ~190 parameter files means any rename of a namelist parameter in `src/` requires a grep over all `problems/**/problem.par*`; there is no automatic check that each parameter file parses (only the ones run by CI do).
