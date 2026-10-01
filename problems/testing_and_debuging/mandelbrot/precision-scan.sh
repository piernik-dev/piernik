#!/bin/bash
set -euo pipefail

# Run one polar Mandelbrot setup sequentially at increasing FP precision.
# Best practiced on problem.par.polar.straight_scepter-quad_precision_test

run_scan() {
	local run_id="$1"
	local precision="$2"

	./piernik -n "&OUTPUT_CONTROL problem_name='mandelbrot_scan', run_id='${run_id}' / &PROBLEM_CONTROL precision='${precision}' /" \
		| tee "out.${run_id}"
}

# Numeric run IDs sort in the same order as precision. Quad complex uses quad
# precision too; it is included last to compare the alternate implementation.
run_scan 001 float
run_scan 002 double-float
run_scan 003 double
run_scan 004 extended
run_scan 005 double-double
run_scan 006 quad
run_scan 007 "quad complex"

../../visual/pvf.py mandelbrot_scan_???_0000.h5 -d mand -D gist_ncar -z 1,5 -m 2
