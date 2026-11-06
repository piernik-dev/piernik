# src/IO – log/tsl output, HDF5 dumps and restarts

Parent: [../AGENTS.md](../AGENTS.md).

## Top level
- `dataio_pub.F90` – **the most widely used module of the code**: `die`, `warn`, `printinfo`, `msg`, verbosity levels (`V_LOG`, ...),
  `code_progress`, namelist handler (`nh`), output parameters and the run id. Almost every module `use`s it.
- `dataio.F90` – (2200 lines) output scheduling and control: namelists `END_CONTROL`, `RESTART_CONTROL`, `OUTPUT_CONTROL`;
  `.log` and `.tsl` (time series) files, `write_data`, `check_log/check_tsl`, user messages (`user_msg_handler`, e.g. stop/dump requests via file),
  walltime-based termination.
- `dataio_user.F90` – hooks for user-supplied additions to tsl/log/HDF5 output (set from `initproblem.F90`).

## hdf5/ (everything pulled by `HDF5`, the default; `I_KNOW_WHAT_I_AM_DOING` disables it)
- Low level: `helpers_hdf5` (overloaded create-attribute/dataset wrappers), `set_get_attributes`, `read_attribute`, `common_hdf5` (shared by dumps/restarts/slices).
- Outputs: `data_hdf5` (single-precision data dumps, GDF-compatible via `grid_data_format` = module `gdf`, readable by yt),
  `restart_hdf5` (front-end) → `restart_hdf5_v2` (current, AMR-capable) and `restart_hdf5_v1` (legacy reader and writer; writes when `use_v2_io = .false.`).
- `particles_io` (`NBODY && HDF5`), `io_debug` (`DEBUG`; namelist `PIERNIK_IO_DEBUG`).
- Parallel-HDF5 collectives: every rank must take part in every dataset call, or the run deadlocks.

## Gotchas
- Changing the layout of restart files requires a version bump and a reader for the old one (see `python/output_upgrade.py`, `bin/fix_restart_*`).
- New field in dumps: register a name in the `*_hdf5` module's variable list and keep names consistent with the Python tools in `python/`.
- Message prefix rule applies: `[dataio:...]`, `[common_hdf5:...]` – the file name, see `bin/checkmessages.sh`.
