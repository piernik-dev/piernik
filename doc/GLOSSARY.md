# Glossary

Short vocabulary of Piernik; reading `AGENTS.md` files and the code is much easier with it. Where a term is defined: file in `src/`.
Terms marked (?) are my best understanding – confirm/correct when you know better.

## Grid and data
| Term | Meaning |
|------|---------|
| `cg` | "computational grid" block: a rectangular chunk of cells plus guardcells, owned by one MPI rank (`grid/grid_container.F90`, type `grid_container`) |
| `cg_list`, `cg_level`, `leaves`, `all_cg` | lists of `cg`s; `leaves` = the finest blocks covering the domain (standard iteration target), `all_cg` = all blocks on all levels (`grid/cg_list*.F90`, `cg_leaves`, `cg_list_global`) |
| `nb` | number of guardcells (boundary cells) around the physical domain; same in all directions (`grid/domain.F90`) |
| `dom` | the global domain description (sizes, boundary types, decomposition) |
| `is, ie, js, je, ks, ke` | index ranges of the physical (non-guardcell) cells of a `cg`; the `RNG` macro in `base/piernik.h` expands to them |
| `qna`, `wna` | registries of *named arrays* – 3D (`q`) and 4D (`w`) fields stored in `cg%q(i)%arr`, `cg%w(i)%arr`; `ind(name)` returns the index (`base/named_array*.F90`) |
| `cg%u`, `cg%b` | fluid state (conserved variables) and magnetic field, stored as `w` arrays |
| `flind` | table of indices of fluid components inside `cg%u` (`fluids/fluidindex.F90`) |
| AMR / refinement | adaptive mesh refinement: blocks on several levels; criteria in `grid/refinement/` (URC = user/unified refinement criteria (?)) |
| prolongation / restriction | interpolating data from a coarse to a finer level / averaging fine to coarse; the order is set by `ord_prolong` (`O_INJ`, `O_I2`, `O_I3`, ... in `constants.F90`) |
| vital array | field that must be kept consistent when blocks are refined/moved (flag in `reg_var`) |
| guardcell exchange / `all_bnd` | filling guardcells from neighbours (MPI) and external boundaries (`base/all_boundaries.F90`) |
| `BND_*` | boundary types (`BND_PER`, `BND_REF`, `BND_OUT`, ...) chosen with `bnd_xl ...` in `BASE_DOMAIN` |
| load balancing | redistributing blocks among ranks (`grid/load_balance*`, `cg_cost*`) |

## Numerics
| Term | Meaning |
|------|---------|
| RTVD | relaxing total-variation-diminishing scheme; oldest solver, dimensionally split (`scheme/rtvd_split`) |
| Riemann solver | MUSCL-type scheme using approximate Riemann solvers (HLLC, HLLD), split or unsplit (`scheme/Riemann_Solver`, `unsplit_schemes`); selected via `which_solver` |
| split / unsplit | 1D sweeps in turn / multi-dimensional update in one step |
| CT | constrained transport – keeps ∇·B = 0 using face-centred B (`scheme/ct.F90`) |
| HDC | hyperbolic divergence cleaning – the alternative ∇·B treatment for Riemann solvers |
| CFL | Courant–Friedrichs–Lewy number: time-step limit (`base/timestep*`) |
| multigrid | solver for the Poisson equation (self-gravity) and CR diffusion (`multigrid/`, `gravity/multigrid_*`) |
| multipole | multipole expansion of the mass distribution, giving external boundary conditions for the Poisson solver (`gravity/multigridmultipole.F90`) |
| FARGO | orbital-advection scheme for rotating discs (`fargo/`) |
| shearing box | local model of a rotating disc; boundary conditions in `shear/` |

## Physics components (preprocessor symbols from `piernik.def`)
`IONIZED` (ionised gas, implies `MAGNETIC`), `NEUTRAL` (neutral gas), `DUST` (pressureless dust fluid), `COSM_RAYS` (CR fluid in the diffusion–advection approximation),
`CRESP` (spectral CR model: momentum bins per cell), `GRAV`/`SELF_GRAV`/`MULTIGRID` (gravity), `NBODY` (particles), `SHEAR`, `RESISTIVE`, `SN_SRC` (supernovae).
Authoritative list and implications: `doc/general/precompiler.F90` (may be incomplete) and `src/base/piernik.h`.

## Run, build and test
| Term | Meaning |
|------|---------|
| problem | directory in `problems/` = initial conditions + parameters (see [../problems/AGENTS.md](../problems/AGENTS.md)) |
| `piernik.def` | compile-time choices for one problem |
| `problem.par` | run-time namelists (`BASE_DOMAIN`, `AMR`, `OUTPUT_CONTROL`, ...) |
| `obj*/`, `runs/` | generated build and run directories |
| `.setuprc` | local default options for `./setup` |
| `pulled by` | comment in a source file listing the symbols that make `./setup` compile it |
| QA | `make QA`: licence, message-format, style and Fortran-convention checks |
| artifacts | `make artifacts`: build/run many configurations and check basic results |
| gold test | regression test comparing a run with a stored reference from a given commit (`jenkins/`) |
| CI | `make CI -j` = QA + artifacts + gold |
| tsl / log / res / h5 | time-series log, text log, restart file, HDF5 snapshot written by `dataio` (`IO/`) |
| PPP | precise parallel profiling – built-in tracing of events (`base/ppprofiling.F90`) |
| hook | procedure pointer through which a problem customises the code (`base/user_hooks.F90`, `IO/dataio_user.F90`) |
| units | code units chosen by `UNITS` namelist; `./setup --units` lists available systems (`base/units.F90`) |
