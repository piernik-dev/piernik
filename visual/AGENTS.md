# visual/ – post-processing and visualisation of Piernik output

Parent: [../AGENTS.md](../AGENTS.md). Standalone Python scripts (no CI). Dependencies: `h5py`, `numpy`, `matplotlib`; `pvf`/`h52paraview` also use Paraview/VTK.
Input is the GDF-compatible `*.h5` dump written by `src/IO/hdf5/data_hdf5.F90`/`grid_data_format.F90` or `.tsl` time-series files.

- `plot_tsl.py` – plot columns of a `.tsl` file against each other/time (`-f`, `-x`, `--xlim`, `--ylim`); python 3-compatible.
- `pvf.py` (+ `pvf_settings.py`, `plot_compose.py`, `plot_utils.py`, `read_dataset.py`) – plots slices/projections of dumps (options incl. scattering particles onto slices);
  settings in `pvf_settings.py`. `pvf.py` is python3; `plot_compose`/`pvf_settings` use `#!/usr/bin/python`.
- `h52paraview.py` – converts dumps (also AMR) to `.vthb` for Paraview; options `-i/-o/-j/-v`.
- Other tools reading Piernik data: `../python/read_h5.py`, `../python/crs_*.py` (CRESP), `yt` can read the files directly (GDF).

## Observations
- Scripts import each other by bare module name, so they must be run from this directory or with `PYTHONPATH=visual`.
- No sample data/tests; a tiny `sedov` dump committed under a test directory would allow smoke tests of the scripts.
- File format changes (names of datasets/attributes in `IO/hdf5`) silently break these scripts; grep here when changing them.
