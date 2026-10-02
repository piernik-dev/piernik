# Mandelbrot test problem

This problem generates Mandelbrot escape-time fields and exercises AMR refinement. It supports Cartesian and logarithmic-polar views, plus scalar, complex, double-float, and double-double precision modes. If `precision` is omitted, the calculation defaults to `double`.

## Main setups

| Parameter file | View and precision | Iteration limit | Notes |
| --- | --- | ---: | --- |
| `problem.par` | Cartesian, double | 1,000 | General overview of the set. |
| `problem.par.polar.straight_scepter-quad_precision_test` | Log-polar, quad | 300 | Deep zoom used by `precision-scan.sh`. |
| `problem.par.polar.Feigenbaum-Myrberg` | Log-polar, extended | 100,000,000 | Very expensive; intended to resolve the high-period structure near the Feigenbaum-Myrberg point. |

## Extra presets

The `extras/` directory contains additional viewports. Precision below is the configured value; presets that omit `precision` use the default `double`.

| Parameter file | View | Precision | Iteration limit |
| --- | --- | --- | ---: |
| `problem.par.zoom_1` | Cartesian zoom | double | 1,000 |
| `problem.par.zoom_2` | Cartesian zoom | double | 10,000 |
| `problem.par.zoom_3` | Cartesian zoom | double | 1,000 |
| `problem.par.zoom_4` | Cartesian zoom | double | 1,000 |
| `problem.par.cartesian.6a6` | Deep Cartesian zoom | quad | 50,000 |
| `problem.par.polar.1-2-3-4-5-6` | Log-polar exploration | double | 10,000,000 |
| `problem.par.polar.2_4_8_16_32` | Log-polar exploration | quad | 200,000 |
| `problem.par.polar.3+3-3+3-3+3` | Log-polar exploration | extended | 3,000,000 |
| `problem.par.polar.3-3-3-3-3-3` | Log-polar exploration | double | 300,000 |
| `problem.par.polar.3a4a5a6a7a8a9` | Log-polar exploration | quad | 500,000 |
| `problem.par.polar.3x3x3x3x3x3` | Log-polar exploration | double | 300,000 |
| `problem.par.polar.6a6` | Log-polar exploration | quad | 50,000 |
| `problem.par.polar.wide_path` | Log-polar exploration | quad | 1,000 |

The preset filename identifies the view; inspect its `xmin`/`xmax`, `ymin`/`ymax`, and center coordinates for the exact region. High iteration limits can make runs substantially slower.

## Output fields

- `mand`: log-like escape-time measure; with smooth coloring enabled, escaped points receive a continuous value rather than the raw iteration count.
- `real`, `imag`: real and imaginary components of the final iterate.
- `dist`: logarithm of the magnitude of the final iterate.
- `ang`: argument of the final iterate, in radians.
- `level`: AMR refinement level.

The `c_polar` parameter biases the displayed escape-time measure by the horizontal coordinate; it does not change the Mandelbrot recurrence.

## Precision scan

Set up the run with `problem.par.polar.straight_scepter-quad_precision_test` as its `problem.par`, then run `precision-scan.sh` from the generated run directory containing `./piernik`. The script executes the modes sequentially with sortable IDs `001` through `007` and uses `../../visual/pvf.py`, relative to that run directory, to visualize the resulting `mandelbrot_scan_???_0000.h5` files. Quad complex is included after quad to compare the alternate implementation at the same nominal precision.
