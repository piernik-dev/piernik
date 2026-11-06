# Cleanup and improvement backlog

Collected while writing the `AGENTS.md` files (see [../AGENTS.md](../AGENTS.md)); details sit in the per-directory files.
Items are grouped by **difficulty** and by **blast radius** (what else can break). Nothing here has been done unless noted.
Every item should be a separate, reviewable pull request. Always finish with `make QA`, and for anything above level A also `make CI -j`.

Legend – blast radius: **none** (docs/ignored files), **local** (one script/tool, not part of the build), **build** (can break `./setup` or compilation for some users/compilers),
**results** (can change numerical results/gold tests), **history** (file moves: `git blame`, open PRs, docs, user scripts).

## A. Easy, no or local side effects (good first tasks, also for AI-assisted newcomers)
| # | Task | Radius | Notes |
|---|------|--------|-------|
| A1 | Remove or port orphaned/Python 2 scripts: `make_test.py`, `python/{plot,drive_via_spawn,repo_crawler}.py`, `bin/{generate_public.py,validate.py,resdiff,fix_restart_1.17_to_1.18.py}` | local | check with `grep -r <name>` first; `repo_crawler` needs private SVN |
| A2 | Add `python3 -m py_compile` of all tracked scripts to `make QA` (guard against silent rot) | local (QA) | do after A1, otherwise QA fails |
| A3 | Unify shebangs to `#!/usr/bin/env python3` | local | macOS has no `python` |
| A4 | Fix `Makefile` header (still mentions `pep8` as "obsoleted"), update `PROJECT_NUMBER` in `piernik.doxy` | none | |
| A5 | `ChangeLog` (frozen 2011) / `NEWS` / `README`: replace with real pointers (GitHub releases, wiki) or delete | none | `COPYING`/`AUTHORS` stay |
| A6 | Delete `doc/Doxyfile.stable` (unused); decide on `GENERATE_LATEX = YES` | none | |
| A7 | Doc hygiene in `src`: missing `\brief` (`piernik.F90`, `fluids_pub`, `cresp_variables`, `fluxlimiters`, `unsplit_*`), wrong brief in `base/randomization.F90` | none | comments only |
| A8 | Fill `info` for problems that lack it; fix `testing_and_debuging` spelling only together with B-level rename (see C3) | none | `./setup --problems` |
| A9 | Remove personal/obsolete compiler configs (`toy`, `gnu47`, `lothlorien`, `varken*`) – ask maintainers first | local | someone may still use them; they are only chosen with `-c` |
| A10 | Bring `doc/general/precompiler.F90` up to date (derive symbol list with grep) and generate a namelist list | none | |
| A11 | Check that `python/maclaurin.py` still works with the current I/O (header claims v1.x) and update the header | local | used by `make Maclaurin` |

## B. Moderate: self-contained but touching shared tooling
| # | Task | Radius | Notes |
|---|------|--------|-------|
| B1 | Port or remove problem helper scripts that use `yt.mods` / `from pylab import *` (8 files, e.g. `mcrtest/piernik_problem.py`) | local | needs yt for testing |
| B2 | Make `bin/bash_completion.sh` generate its option list from `./setup --help` | local | |
| B3 | Deduplicate `compilers/*.in` with a shared `common.in` (`default.in` = `h5pfc_gnu.in`) | build | risk: different flag order on some machines; test with the CI configs (`jenkins.in`, `default.in`) and a `diff` of generated `Makefile`s |
| B4 | Make `compilers/default.in` Apple-aware for Homebrew on Apple silicon (`/opt/homebrew`) | build | needs a Mac tester; see the macOS notes in `python/AGENTS.md` |
| B5 | Re-check `setup` portability: `gcc -dM -E -` (OS detection), `pkg-config` (fftw3, fortrancl) on macOS; finish the `cpp` fix after feedback from the Mac user | build | gate: Makefiles for `sedov`, `jeans`, `maclaurin` identical before/after (sorted lines, `PYTHONHASHSEED=0`) |
| B6 | Move inline shell test logic and hard-coded expected numbers out of the root `Makefile` into scripts in `bin/`/`jenkins/` | build (CI) | keep target names (`make CI`, `QA`, `artifacts`, `gold`) stable; compare outputs before/after |
| B7 | Replace `eval $(sed ...)` config parsing in `jenkins/gold_test.sh` with a safer parser; make shallow-clone failure message explicit | local (CI) | |
| B8 | Add a smoke test for `visual/` scripts with a tiny committed HDF5 dump | local | needs a small data file (<100 kB) |
| B9 | Add a "parameter files parse" check: run `./setup <problem> -p <file> -n` and a `mpirun -n 1 ... tend=0` for all `problem.par*` (at least those with `.build`) | none (new test), but long | detects stale parameter names after namelist changes |
| B10 | Regenerate or document golds workflow; remove `bin/regenerate_golds.sh` (private paths) | local | superseded by commit-based gold test |

## C. Hard: wide side effects
| # | Task | Radius | Notes |
|---|------|--------|-------|
| C1 | Reorganise `src/` thematic split (see "Proposed cleanups" in [../src/AGENTS.md](../src/AGENTS.md)): `grid/balance`, `base/{ppp,containers}`, gathering multigrid into `multigrid/`, merging single-file dirs | history | `setup` globs `src/**`, so builds keep working, but **module names must remain unique** and `pulled by` lines must travel with files; batch moves in one PR, update docs and `AGENTS.md` at once, announce before open PRs conflict |
| C2 | Move physics out of `base/` (`timestep*`, `all_boundaries`, `coriolis`, ...) to restore the layering | history + build | cyclic `use` dependencies may appear; `make dep` shows the graph |
| C3 | Rename inconsistent problems (`testing_and_debuging`, `stream*`, `2body/3body` vs `N-body`) | history | breaks `.setup.call` files, gold configs (`jenkins/gold_configs`), scripts and wiki; do with deprecation notes |
| C4 | Fix `pulled by` mismatches (`particle_types`/`particle_func` by `GRAV` vs `NBODY`, `particles_io`) | build | can change which files are compiled for some `piernik.def` combinations – run `make allsetup` and compare the generated file lists before/after |
| C5 | Rewrite `setup` file-selection (`pulled by`, `overrides`, regex on comments) into a declarative, testable form; add unit tests for the parser | build | highest risk to every user; keep Makefile output identical (verify on all problems with `make allsetup`) |
| C6 | Unify the `piernik.def` / `piernik.h` symbol system (many derived symbols, undocumented combinations) and validate unsupported combinations early in `setup` | build + results | |
| C7 | Remove Doxygen's duplication by moving namelist documentation next to declarations and generating pages | none–local | big editorial effort |
| C8 | Retire the legacy I/O v1 remnants (`python/maclaurin.py` assumptions, restart converters) once no data needs them | results | needs maintainers' decision on supported file versions |

## D. Needs a decision first (not tasks yet)
- Which problems are still maintained? Mark unmaintained ones (e.g. in `testing_and_debuging/` and the 8 problems with broken Python helpers) in `info` or move to an `attic/`.
- Should `make CI` also build a macOS/ifx configuration? (currently only gfortran on Linux in `jenkins.in`).
- Which compilers/MPI/HDF5 versions are officially supported? Write them in `INSTALL.md`; this also allows dropping workarounds (`-fallow-argument-mismatch`, `NO_F2018`, `FORBID_F08`).

## Suggested order
A1 → A2 → A4–A7 (one afternoon, almost no risk) → B5 (waiting for Mac feedback) → B9 (safety net) → B3/B6 → then C4 → C1/C3 together, once B9 protects against breakage.
