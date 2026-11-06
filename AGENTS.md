# Piernik – context for developers and AI agents

Piernik is a parallel (MPI + HDF5) Fortran 2008 MHD code with AMR, multigrid self-gravity,
cosmic rays (incl. the CRESP spectral model), dust/neutral/ionized fluids and N-body particles.
Documentation entry points: `README`, `INSTALL.md`, `doc/general/*.F90` (Doxygen pages: namelists,
precompiler symbols, how to write `initproblem.F90`, `new_source_file.F90` template).

## How the context files are organised (keep this scheme when adding more)
Context is **hierarchical**: each file is short, describes only its own directory and points to
neighbours. Tools that support nested `AGENTS.md` load the closest one, so an agent working in
`src/grid` reads this file + `src/AGENTS.md` + `src/grid/AGENTS.md` and nothing else.
Do not duplicate information between levels; link instead.

| File | Describes |
|------|-----------|
| `AGENTS.md` (this) | project-wide facts, build/test, conventions |
| [python/AGENTS.md](python/AGENTS.md) | the `./setup` script (`python/piernik_setup.py`) and the other Python tools |
| [bin/](bin/AGENTS.md), [compilers/](compilers/AGENTS.md), [jenkins/](jenkins/AGENTS.md), [visual/](visual/AGENTS.md), [benchmarking/](benchmarking/AGENTS.md) | helper scripts, build configs, CI/gold tests, post-processing, benchmarks |
| [src/AGENTS.md](src/AGENTS.md) | source tree map, layering, preprocessor symbols, small directories |
| [problems/AGENTS.md](problems/AGENTS.md), [doc/AGENTS.md](doc/AGENTS.md) | problem anatomy and test problems; Doxygen pages and docs |
| `src/{base,grid,fluids,scheme,IO,gravity,multigrid,magnetic,particles}/AGENTS.md` | thematic areas |

Guideline for new files: one per directory that has ≥ ~8 modules or a non-obvious design; keep each
below ~100 lines; describe *roles and relations*, not every routine (the code has Doxygen `\brief`s).
Run `grep -rn "\\brief" <file>` or `./setup --problems` rather than copying lists that rot quickly.
When you move/rename modules, update the corresponding file in the same commit.

## Build and run
- `./setup <problem> [-c compiler] [-d SYMBOLS] [-n]` creates `obj*/` and `runs/<problem>/`; see
  [python/AGENTS.md](python/AGENTS.md). `./setup --problems` lists problems (`problems/*/info`).
- Problem = directory in `problems/` with `piernik.def` (preprocessor symbols selecting physics),
  `initproblem.F90` (ICs, user hooks), `problem.par` (namelists, run-time parameters), optional `info`.
- Which source files get compiled is decided by `setup` from `! pulled by SYMBOL` comments plus
  the symbols in `piernik.def` / `src/base/piernik.h`. Run-time choices (solver, BCs, ...) are
  namelists in `problem.par`.
- Compiler configs: `compilers/*.in`. Typical dependencies: gfortran/ifx, MPI, parallel HDF5 (`h5pfc`), FFTW.

## Tests
- Formal suite: `make CI -j` = `QA` (license headers, error messages, pycodestyle, qa.py) +
  `artifacts` (setup+build of many configurations) + `gold` (regression against stored results).
  See the top of `Makefile` and `make help`. Run groups separately: `make QA`, `make artifacts -j`, `make gold -j`.
- Standard small problems: `sedov`, `maclaurin`, `jeans`, `otvortex`, `crtest`, `mcrtest`, `3body_orig`,
  `streaming_global`, `advection_test`. Quick sanity check after changing the code:
  `./setup sedov -n` (and a build), then `make QA`.
- `jenkins/` holds the gold-test scripts and configs used by the Makefile.

## Conventions
- Every module should have a license header (checked by `bin/check_license_headers.sh`) and a Doxygen
  `\brief`; start new files from `doc/general/new_source_file.F90`. The modules which lack Doxygen comments, should be updated to include them.
- Messages: `die`, `warn`, `printinfo` from `dataio_pub`, text prefixed with `"[file_basename:routine]"`
  (enforced by `bin/checkmessages.sh`; this is the *file* name, which may differ from the module name).
- `bin/qa.py` (part of `make QA`) enforces: `implicit none`, explicit `public`/`private`, `use ... only:`,
  `::` in declarations, no tabs, no implicit `save` through initialised locals, few magic integers;
  suppress a warning only with a justified `QA_WARN` comment.
- MPI: use the wrappers in `src/base/mpi` (`piernik_MPI_*`) instead of raw `MPI_*` calls where one exists.
- Preprocessor symbols belong to `piernik.def`; derived symbols are set in `src/base/piernik.h`.
- Python scripts must pass `pycodestyle` (long lines ignored).
- Generated artefacts (`obj*/`, `runs/`, `*.png`, `tmp*.f90`) should not be committed.

## Known technical debt (details in the per-directory files)
- Python 2 scripts that no longer run (`python/{plot,drive_via_spawn,repo_crawler}.py`, `bin/{generate_public,validate,resdiff}`, ...) and are not covered by CI.
- Compiler configs are copy-pasted (`default.in` = `h5pfc_gnu.in`), several are personal/historic (`toy`, `gnu47`, `lothlorien`, `varken*`).
- `bin/bash_completion.sh` duplicates the `setup` option list; `bin/regenerate_golds.sh` has private absolute paths.
- Directory split in `src/` (see "Proposed cleanups" in [src/AGENTS.md](src/AGENTS.md)).
- Some tests/documentation assume Linux paths (`/usr/include`, Fedora module dirs); macOS (Homebrew in `/opt/homebrew`) needs care.
- Root-level files: `ChangeLog` frozen in 2011 (SVN era), `NEWS` only points to a wiki, `README` is a stub; `make_test.py` is an
  orphaned Python 2 script not referenced by Makefile/CI (delete or port); `piernik.doxy` is a 2.5 k-line default dump.
- `Makefile` is large and contains test logic (inline shell, hard-coded expected numbers) – candidates for moving into `bin/`/`jenkins/` scripts.
  Its header still advertises `pep8` as obsolete; `ALLOBJ = $(wildcard obj*)` and the `*:` catch-all rule make stray `obj*` dirs affect `make clean`/`all`.
- `problems/` and `doc/`: see [problems/AGENTS.md](problems/AGENTS.md) and [doc/AGENTS.md](doc/AGENTS.md) (old yt scripts, stale docs, inconsistent names).

## Root files at a glance
`Makefile` (CI, QA, local helpers: `resetup`, `allsetup`, `ctags`, `dep`/`view_dep`, `doxy`, modifiers `RS=1`/`CL=1`), `INSTALL.md` (install guide, best human-oriented overview),
`.gitattributes` (`merge=ours` for `.setuprc`, `runs/check.sh`), `.setuprc[.$HOSTNAME]` (per-machine default setup options, git-ignored),
`setup` (symlink to `python/piernik_setup.py`), `AUTHORS`, `COPYING` (GPL/licence text; headers checked by `make QA`).
