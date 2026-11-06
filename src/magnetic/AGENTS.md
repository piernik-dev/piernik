# src/magnetic – magnetic field tools

Parent: [../AGENTS.md](../AGENTS.md). `MAGNETIC` is implied by `IONIZED` unless `NONMAGNETIC` is defined.

- `divB.F90` (module `div_B`) – ∇·B diagnostics (`print_divB_norm`) and cleaning helpers.
- `bfc_bcc.F90` – face-centred ↔ cell-centred B conversion (B is stored face-centred for CT, `cg%b`).
- `magboundaries.F90` – BCs for B, vector potential and EMFs.
- `resist/resistivity.F90` (`RESISTIVE`, namelist `RESISTIVITY`) – resistive dissipation, timestep limit.
- Related but located elsewhere: `scheme/ct.F90` (constrained transport), `scheme/Riemann_Solver/hdc.F90` (divergence cleaning), `scheme/Riemann_Solver/hlld.F90` (MHD Riemann solver).

Tests: `otvortex`, `field_loop_mhd`, `cpaw`, `magnetic_rotor`; gold test `gold_resist`.
