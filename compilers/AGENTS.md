# compilers/ – build configurations

Parent: [../AGENTS.md](../AGENTS.md). Selected with `./setup <problem> -c <name>` (name without `.in`; default `default.in`).
The chosen file is **pasted at the top of the generated `obj*/Makefile`** (see [../python/AGENTS.md](../python/AGENTS.md)),
followed by the templates that run feature tests from `tests/`. It is therefore plain GNU make syntax.

## Variables a config provides
`PROG` (executable name), `F90` (compiler wrapper, normally `h5pfc` or `mpif90`), `F90FLAGS`, `LDFLAGS`, optionally `LIBS`, `CC`/`CFLAGS`
(for `defines.c`), `MAKEFLAGS`, `PRECOMP`. `PIERNIK_DEBUG=1` (setup `--debug`) is available as an `ifeq ($(PIERNIK_DEBUG), 1)` switch.
Mandatory flags for gfortran: `-fdefault-real-8 -ffree-form -fimplicit-none -ffree-line-length-none` (all reals are double via the compiler flag, not `kind=`).

## Files
| File | Purpose |
|------|---------|
| `default.in` (= `h5pfc_gnu.in`, identical copy) | recommended for typical Linux/macOS installs, follows `INSTALL.md` |
| `camk.in`, `camk_paranoid.in`, `varken_camk.in` | developer (Toruń) configs: gfortran, `camk_paranoid` adds very strict checking |
| `ifx.in` | Intel oneAPI (`-DNO_F2018`); for compatibility checks |
| `jenkins.in` | config for the Jenkins server (Fedora, `-O3`, `-j4`); `make CI`/gold tests do **not** pass `-c`, so they use `default.in` unless a `.setuprc` says otherwise |
| `lothlorien.in`, `varken*.in`, `toy.in`, `gnu47.in` | machine- or person-specific and historic (gcc 4.7 / 5.1); candidates for removal |
| `newcompiler` | 2008 perl helper listing configs; obsolete (refers to `../obj`) |
| `tests/` | tiny Fortran programs compiled by the generated Makefile: `mpi_f08.F90`, `mpi.F90`, `F2018.F90`, `mpi_allgatherv_bug.F90` (run!), `mpi_fat_ring_test.F90`. Results become `-DMPIF08`/`-DNO_MPIF08_AVAILABLE`, `-DNO_ALL_MPI_FUNCTIONS_AVAILABLE`, `-DNO_F2018`, `-DFORBID_F08` |

## Notes
- Add a new config by copying `default.in`; do not edit `default.in` for a local setup. Local files are not ignored by git (`*.in`), so do not commit personal configs by accident.
- `default.in` contains Apple-specific branching (`TERM_PROGRAM`, `/usr/local/include`); Homebrew on Apple silicon uses `/opt/homebrew` – a likely source of macOS problems
  (include path, `cpp`; see the macOS note in `python/AGENTS.md`).
- Many warnings from MPI (`-fallow-argument-mismatch`, `-frecursive`) are described in `INSTALL.md`; do not add `-frecursive` (breaks CRESP initialisation with gfortran 14).

## Observations for cleanup
- `default.in` and `h5pfc_gnu.in` are duplicates; `camk.in`/`lothlorien.in`/`jenkins.in` share most lines. A shared `common.in` (included via `include`) would remove the copy-paste.
- Version-named configs (`gnu47`, `toy` with gcc-5.1) are unmaintained.
