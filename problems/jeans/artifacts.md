# Interpreting the results in `artifacts` test

## The default, 64×64, 2D test

The results of the fit depend on (in)accuracy od the solver. Both amplitude and period may be affected.

As of late 2025 (merged PR #598) the values are:

solver            amplitude  (rel. error)  period   (rel. error)  damping
----------------- ---------- ------------- -------- ------------- --------
analytical        3.914e-7    —            13.736    —             —
HLLC              4.175e-7    0.063        14.181    0.031        -0.014463
RTVD              3.919e-7    0.0013       13.719   -0.0012       -0.0063501
Riemann split     3.913e-7   -0.00025      13.815    0.0057        0.011448
Riemann unsplit   3.9129e-7  -0.00027      13.815    0.0057        0.011514

* It seems that the HLLC solver has some major issues. It is not surprising as it wasn't maintained for years and is not planned to be used in production runs.
* Both split and unsplit Riemann solvers behave very similar to each other.
* The Riemann solvers are very accurate at recovering the amplitude while the RTVD is slightly better at oscillation period and damping factor.

## Increased resolution to 256×256

The results are as follows:

solver            amplitude  (rel. error)  period   (rel. error)  damping
----------------- ---------- ------------- -------- ------------- --------
analytical        3.914e-7    —            13.736    —             —
HLLC              4.175e-7    0.063        14.186    0.032        -0.00023135
RTVD              3.9139e-7  -0.000024     13.735   -0.000069      0.01427
Riemann split     3.914e-7   -0.0000085    13.741    0.0003        0.014478
Riemann unsplit   3.9141e-7   0.000018     13.74     0.0003        0.010059

* The errors in HLLC barely changed, which strongly suggest that the implementation suffers from some mistake.
* The amplitude and period errors in RTVD and Riemann solvers improved greatly.
* The damping factor remained roughly the same which may indicate limited temporal accuracy of the self-gravity solver and its coupling via the source terms.
* The Riemann unsplit solver required reduction of `cfl` factor to 0.5, otherwise it failed badly.
