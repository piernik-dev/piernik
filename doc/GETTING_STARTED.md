# Getting started with Piernik

[INSTALL.md](../INSTALL.md) gets the code compiled. This page is the next step: first run, what is where, and how to start changing things.
Developers who want to use AI assistants: see [AI_DEVELOPMENT.md](AI_DEVELOPMENT.md).

## 1. First run (15 minutes)
```
./setup sedov                 # builds obj/ and prepares runs/sedov/
cd runs/sedov
mpirun -n 4 ./piernik         # tend = 1e-3 in problem.par, so it is short
```
Results appear in the run directory: HDF5 snapshots (`*.h5`, every `dt_hdf`), a text log (`*.log`), and time-series (`*.tsl`). Restart files (`*.res`) are written only when restart output is enabled with a positive `dt_res`.
`./setup --problems` lists all problems with one-line descriptions; `./setup --help` lists options.

## 2. The three levels of configuration
| What | Where | When it takes effect |
|------|-------|----------------------|
| Physics modules, solver variants (`IONIZED`, `GRAV`, `COSM_RAYS`, ...) | `problems/<name>/piernik.def` and `./setup -d SYMBOL` | compile time (re-run `./setup`) |
| Domain size, resolution, BCs, `tend`, output intervals, solver options | `problems/<name>/problem.par` (copied to `runs/<name>/`) | run time: edit the copy in `runs/`, no recompilation |
| Initial conditions and user hooks | `problems/<name>/initproblem.F90` | compile time |

Edit `runs/<name>/problem.par` freely – `./setup` overwrites it only when you call it again (use `-k` to keep your copy).
Parameter meanings: Doxygen pages (`make doxy`, then `doc/html`) and comments next to each `namelist` in `src/`
(`grep -rn "namelist /" src`).

## 3. Everyday setup options
`-c <config>` compiler config from `compilers/`; `--debug` checks and debug flags; `-n` only prepare the build (cheap sanity check);
`-o <suffix>` separate `obj_<suffix>/` and `runs/<name>_<suffix>/` (keep several variants side by side); `-p problem.par.<variant>` alternative parameter file;
`--last` repeat the last call. Put your habitual options in `.setuprc` (per-machine: `.setuprc.$HOSTNAME`), e.g. `-c mycluster --laconic --linkexe`.
Source `bin/bash_completion.sh` for tab completion of options and problem names.

## 4. Typical tasks
| I want to... | Do this |
|--------------|---------|
| change resolution / end time / output frequency | edit `runs/<name>/problem.par` (`BASE_DOMAIN`, `END_CONTROL`, `OUTPUT_CONTROL`) |
| run my own initial conditions | copy the closest problem dir, edit `initproblem.F90` and `problem.par`; see [problems/AGENTS.md](../problems/AGENTS.md) and `doc/general/initproblem.F90` |
| add/remove physics | change `piernik.def` of your problem; see symbols in `doc/general/precompiler.F90` and `src/base/piernik.h` |
| restart a run | rerun `./piernik` in the run directory; it picks up the newest `*.res` (see `INSTALL.md`/`problem.par` for `restart` options) |
| look at results | `visual/` scripts, yt, or VisIt/ParaView on the `*.h5` files (see [visual/AGENTS.md](../visual/AGENTS.md)) |
| find a module/routine | [src/AGENTS.md](../src/AGENTS.md) has a map of the sources and the program flow |
| debug | `./setup <problem> --debug -o dbg`; run under `mpirun -n 1` first |

## 5. Contributing
1. Fork on GitHub, work on a branch (never commit directly to `master`).
2. Before a pull request run `make QA` (fast: licence headers, message format, style, Fortran conventions) and, for anything touching the numerics, `make CI -j 3`.
3. Conventions (messages, MPI wrappers, `only:` lists, licence header): [AGENTS.md](../AGENTS.md#conventions) – they are checked automatically.
4. Changing the results of a standard test deliberately? Say so in the PR; gold tests will differ (see [jenkins/AGENTS.md](../jenkins/AGENTS.md)).

## 6. Common pitfalls
- Missing/old MPI or HDF5 wrappers: `./setup` fails in the compiler test – see the "Other systems" and Intel sections in `INSTALL.md`.
- macOS: Homebrew lives in `/opt/homebrew`; compilers' config `default.in` may need paths adjusted; `setup` needs a working `cpp` (any, but not an empty `$CPP`).
- Module names must be globally unique (all sources are flattened into one build directory).
- A problem that compiled yesterday and not today: `make allsetup` (sets up every problem) quickly shows what your change broke.
