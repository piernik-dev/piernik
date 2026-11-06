# Piernik coding instructions

Keep these project-wide instructions concise. Before changing code, read the root
[AGENTS.md](../AGENTS.md) and the closest relevant directory-level `AGENTS.md`; those files
are the source of truth for detailed architecture, conventions, and subsystem-specific tests.
Do not duplicate their inventories here.

## Project model

- Piernik is a Fortran 2008 simulation using MPI, parallel HDF5, and AMR. Numerical and
  distributed behavior must be preserved unless the requested change explicitly alters it.
- `./setup <problem>` assembles and builds an independent `obj/` directory from `src/` and
  the selected problem. Source selection depends on `piernik.def`, derived symbols in
  `src/base/piernik.h`, and source comments such as `! pulled by SYMBOL`; module dependencies
  are discovered from `use` and `#include`.
- `piernik.def` selects compile-time features; `problem.par` supplies run-time namelist
  parameters. Do not substitute one mechanism for the other.
- The setup flattens selected modules into one build directory: module names must be globally
  unique, regardless of their source subdirectory.

## Implementation rules

- Follow the Fortran checks in `make QA`: use `implicit none`, explicit `public`/`private`,
  `use ... only:`, `::` declarations, and avoid tabs, implicit `save` from initialized locals,
  and unexplained magic integers.
- New Fortran modules need the standard license header from
  `doc/general/new_source_file.F90` and a Doxygen `\brief`. Use `die`, `warn`, and `printinfo`
  from `dataio_pub`; message prefixes identify the **file** (`[file_basename:routine]`).
- Use the `piernik_MPI_*` wrappers in `src/base/mpi` where available. Use `fluidindex` rather
  than hard-coded fluid-state indices.
- Preserve the default-real-8 compiler semantics. Treat changes to numerical algorithms,
  initial conditions, defaults, MPI collectives, AMR, or restart formats as behavior changes:
  justify them and run the relevant regression tests. Restart-layout changes need backward
  compatibility and a versioned reader.
- Keep changes scoped. Do not commit generated build/run output (`obj*/`, `runs/`, temporary
  Fortran files, or generated images); do not remove or overwrite user run data.

## Build and verification

- For a quick setup check, run `./setup <problem> -n`; for a build, run `./setup <problem>`.
  The setup command must be run from the repository root.
- Run `make QA` for Fortran conventions, license/message checks, and Python style. Python
  changes can also be checked with `make pycodestyle`.
- For numerical or build-system changes, build and run a relevant small problem, then run
  the applicable `make artifacts` and/or named gold test. Use `make CI -j` for broader
  validation when practical; report any tests that could not be run.

## Directory guides and references

- [Setup and Python tools](../python/AGENTS.md), [compiler configurations](../compilers/AGENTS.md)
- [Source map and layering](../src/AGENTS.md); each source area has its own `AGENTS.md`
- [Problems and problem tests](../problems/AGENTS.md), [gold tests and CI](../jenkins/AGENTS.md)
- [QA and helper scripts](../bin/AGENTS.md), [documentation](../doc/AGENTS.md)
- [Visualisation tools](../visual/AGENTS.md), [benchmarks](../benchmarking/AGENTS.md)
- [AI development workflow](../doc/AI_DEVELOPMENT.md), [getting started](../doc/GETTING_STARTED.md),
  [installation](../INSTALL.md)
