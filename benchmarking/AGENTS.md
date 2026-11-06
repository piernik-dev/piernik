# benchmarking/ – single-node scaling tests

Parent: [../AGENTS.md](../AGENTS.md). Full description for users: [README.md](README.md); keep this file to the points that matter when editing.

- `piernik_bench.sh` – runs a set of problems (sedov, crtest, maclaurin) once, for 1..N threads; prints a timing table (`-c file` selects `benchmarking/<file>.in` as the compiler config).
- `standard_set_of_tests.sh <results_dir>` – 3 repetitions of the standard size + two enlarged ones.
- `benchmarking.in`, `benchmarking_ifx.in` – compiler configs (`-O3 ... -w`, no debug checks); used via setup's `-c` option, which accepts a path relative to `compilers/`
  (that is why the script uses `../benchmarking/...`).
- `problem.par.{sedov,maclaurin,crtest}` – parameter files with the sizes of the benchmarks; changing them invalidates comparison with old results.
- Plot results with `bench_plot.py` from the external repository `gawrysz/Piernik-benchmarks` (not in this repo).

Notes: Do not run on a loaded/shared machine; the scripts use all hardware threads. `benchmarking.in` is another copy of the gfortran flags (see compilers cleanup notes).
