# doc/ – documentation sources

Parent: [../AGENTS.md](../AGENTS.md). Documentation is Doxygen-based and lives mostly *inside the Fortran sources* (`\brief`, `!<` comments).

- `general/*.F90` – Doxygen "pages" (not compiled): `main.F90` (main page), `namelists.F90` (all namelists), `precompiler.F90` (preprocessor symbols),
  `initproblem.F90` (contract of `initproblem.F90` and user hooks), `new_source_file.F90` (**template of a new source file; its first 26 lines are the licence
  header that `bin/check_license_headers.sh` compares with all sources**).
- `../piernik.doxy` – the Doxygen config actually used by `make doxy` (`INPUT = src doc/general`, output to `doc/html` and `doc/latex`, both git-ignored; `GENERATE_LATEX = YES`).
- `Doxyfile.stable` – older config, unused by the Makefile.
- `images/*.eps` – figures referenced by the Doxygen pages (EPS only; no PNG/SVG for web output). `weno.gnuplot` – plot script for the WENO notes.

## Human-oriented guides (Markdown)
`GETTING_STARTED.md` (first run, configuration levels, contributing), `HOWTO.md` (recipes: add parameter/hook/field/output/BC/...), `GLOSSARY.md` (vocabulary), `AI_DEVELOPMENT.md` (using AI assistants), `BACKLOG.md` (cleanup tasks by difficulty and blast radius). Keep them consistent with `../INSTALL.md` and the root `AGENTS.md`.

## Also read (outside `doc/`)
`../INSTALL.md` (install, MPI/HDF5 pitfalls, Intel compilers, CI usage), `../README` (pointer to the web page), `../Makefile` header (usage of all targets),
`../jenkins/README.md`, `../benchmarking/README.md`, `../problems/jeans/artifacts.md` (interpretation of test results).

## Observations
- Precompiler docs (`precompiler.F90`) list symbols by hand and lag the code: e.g. `NBODY`, `CRESP`, `THERM`, `RANDOMIZE`, `NONMAGNETIC`, `NO_F2018`, `MPIF08`, `FORBID_F08`, `NO_FFT`
  are used in sources/Makefile but only partly documented. Derive the list with `grep -rhoE "#\s*(if|ifdef|ifndef|elif)[^/]*" src | sort | uniq -c`, then update.
- `namelists.F90` (41 lines) is a stub: namelists are documented only next to their declaration in the modules (`grep -rn "namelist /" src`). A generated list would be better.
- `piernik.doxy` is 2.5 k lines of defaults – a minimal config overriding only changed options would be easier to maintain; `PROJECT_NUMBER = 1.0` is stale.
- Doxygen output is not built by CI, so broken `\ref`/`\copydoc` links (e.g. in `initproblem.F90` page) go unnoticed.
- No user manual in the repo; the wiki/website (`piernik.umk.pl`, old `piernik.astri.umk.pl`) is referenced from `README`/`NEWS` and may be outdated.
