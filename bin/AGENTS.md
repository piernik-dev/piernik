# bin/ – developer, QA and debugging helpers

Parent: [../AGENTS.md](../AGENTS.md). Scripts are run from the repository root unless stated otherwise.

## Used by `make QA` (keep them working and fast – they gate every CI run)
- `qa.py` – static checks of the Fortran sources (see the list in the root `AGENTS.md`); also prints warnings with `QA_WARN` escapes. Python 3, needs numpy.
- `checkmessages.sh` – message prefixes `[file_basename:routine]` must match the file they are in.
- `check_license_headers.sh` – first 26 lines of every `.F90` must be identical to `doc/general/new_source_file.F90`
  (a header change means a mass edit of all files; check with this script).
- `implicit_save.sh` – helper for detecting initialised locals (implicit `save`), uses sed on a temporary buffer.
- `generate_public.py` – proposes `public` lists (imports `qa`). **Python 2 syntax**, does not run under python3.

## Testing / regression
- `regenerate_golds.sh` – rebuilds gold files from fixed historical commits; **contains an absolute path to a private disk** (`/raid_hum/...`)
  and old compiler names (`gnu47`): treat as documentation of how golds were produced, not as a portable tool.
  Current gold tests are defined in `../jenkins/gold_configs/` (see [../jenkins/AGENTS.md](../jenkins/AGENTS.md)).
- `gdfdiff`, `gdf_distance`, `resdiff` – compare two HDF5 outputs (GDF dumps / restart). `gdf_distance` is the modern one (py3, argparse, typing);
  `gdfdiff` shebang is `/usr/bin/python` and `resdiff` is python2 (uses `scipy.misc`, now removed).
- `validate.py` – old (python2) validation driver with hard-coded paths; superseded by `make CI`.

## Run-time / debugging
- `gdb_wrap`, `gdb_find_div_in_bt`, `gdb_killall` – catch where a parallel run deadlocked (backtrace per rank, then find the first divergence).
  `gdb_killall` kills by host name only – dangerous on shared machines.
- `extract_problem_par.{sh,py}` – recreate a `problem.par` from `.h5` or `.log` (the `.log` variant also contains command-line overrides).
- `compress_hdf5_3D.sh` – (re)compress 3D dumps (not restarts); `fix_restart_1.17_to_1.18.py` – one-off restart converter (ancient); python2.
- `f90grep`, `f90sed` – grep/sed over all Fortran sources; need `PIERNIK_ROOT` (hard-coded old layout assumptions).
- `bash_completion.sh` – `source` it for completion of `./setup` options and problem names; the option list is **a copy** of the `setup` options
  and must be updated by hand when options change.

## Observations for cleanup
- Python 2 leftovers: `generate_public.py`, `validate.py`, `resdiff`, `fix_restart_1.17_to_1.18.py`, plus `python/{drive_via_spawn,plot,repo_crawler}.py`.
  Port to python3 or delete; they are not covered by CI, so they rot silently. A cheap guard: `python3 -m py_compile` on all scripts in the `pycodestyle`/QA target.
- Mixed shebangs (`python`, `python3`, `/usr/bin/python`, `/usr/bin/env python`) – on macOS/modern Linux plain `python` often does not exist; use `#!/usr/bin/env python3`.
- `bash_completion.sh` option list duplicates `setup`; generate it from `./setup --help` instead.
- No tests for the helpers themselves.
