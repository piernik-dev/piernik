# src/fluids – fluid components, sources, cosmic rays

Parent: [../AGENTS.md](../AGENTS.md). Each component has its own subdirectory (mirrors the preprocessor symbols).

## Top level
- `fluidindex.F90` – **the index map** of the state vector: `flind` (counts, `flind%ion`, `flind%neu`, `flind%dst`, `flind%crs`, ...)
  and the position of each variable in `cg%u(:, :,:,:)`. Everything that loops over fluids uses it.
- `fluidtypes.F90` – compound types describing a fluid (`component_fluid`, `var_numbers`, ...); `fluids_pub.F90` – a few shared switches.
- `initfluids.F90` – initialises all components, tracers and dependent variables (called from `initpiernik`);
  order of calls to `init_*` matters because they register fields and indices.
- `fluxes.F90` – gathers flux components from every fluid for the solvers; `sources.F90` – gathers all source terms
  (gravity, Coriolis, CR, interactions, thermal, ...) applied by the integration schemes.

## Components (pulled by symbol, see `piernik.def`)
| Dir | Symbol | Contents |
|-----|--------|----------|
| `ionized/` | `IONIZED` | `initionized` (also sets up magnetic field when `MAGNETIC`) |
| `neutral/` | `NEUTRAL` | `initneutral` |
| `dust/` | `DUST` | `initdust` |
| `tracer/` | – | `inittracer`, `fluxtracer` (passive scalars) |
| `interactions/` | `FLUID_INTERACTIONS` | drag etc. between fluids, `timestepinteractions` |
| `thermal/` | `THERM` | `thermal.F90` – radiative cooling/heating (implicit energy update) |
| `cosmicrays/` | `COSM_RAYS`, `CRESP` | see below |

Each `init*` reads its namelist (e.g. `COSMIC_RAYS`), sets EOS parameters and registers variables in `flind`.
Isothermal vs adiabatic is a compile-time choice (`ISO`) plus run-time parameters.

## cosmicrays/
- Fluid-like CR energy species: `initcosmicrays`, `cr_data`, `fluxcosmicrays`, `sourcecosmicrays`, `timestepcosmicrays`,
  `crdiffusion` (explicit anisotropic diffusion), `crhelpers`, `multigrid_diffusion` (implicit diffusion; needs `MULTIGRID`).
- **CRESP** (spectral model, `CRESP`): `initcrspectrum` (namelist, spectrum grid), `cresp_crspectrum` (core algorithm, ~2000 lines),
  `cresp_NR_method` (Newton–Raphson solvers for spectrum reconstruction), `cresp_grid` (apply to every cell),
  `timestep_cresp`, `cresp_io`, `cresp_helpers`, `cresp_variables`. Numerically delicate – use the `crtest`/`mcrtest` problems and gold tests
  `gold_CRESP`, `gold_mcrtest`, `gold_mcrwind` after any change. Compiling with `-frecursive` is known to break CRESP initialisation (see `INSTALL.md`).

## Adding a new fluid/field – checklist
1. new symbol in `piernik.def` docs (`doc/general/precompiler.F90`) and a `! pulled by` line in new files;
2. extend `fluidindex` (+ `fluidtypes`), register in `initfluids`;
3. flux in `fluxes`, source in `sources`, dt limit in `base/timestep`;
4. boundaries (`base/fluidboundaries_funcs`, `scheme`), IO names (`IO/hdf5/*`), restart compatibility.
