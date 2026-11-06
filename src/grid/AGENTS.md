# src/grid – domain, grid containers, lists/levels, AMR, load balance

Parent: [../AGENTS.md](../AGENTS.md). 51 modules + `mergeboxes/`, `refinement/`. Largest area of the code; most
algorithms here are MPI-collective – every rank must call them in the same order.

## Domain and decomposition
`domain.F90` (`dom`: global/local box, boundary types, geometry sizes), `decomposition.F90` (cartesian/AMR block decomposition),
`patch_list`, `ordering` (space-filling curve ids), `mergeboxes/` (convert a 3D logical map into few boxes; used for decomposition),
`level_essentials`, `grid_helpers`, `dot.F90` (topology picture for debugging).

## Grid container `cg` (one block of cells + guardcells)
Type chain (each file adds one feature; all abstract except the last):
`grid_container_base` → `_na` (named arrays) → `_bseg` (boundary segments) → `_op` (grad/div/curl/dot) → `_ref` (refinement)
→ `_fcflx` (fine/coarse fluxes) → `_prolong` → `grid_container` (final). Add features in the matching layer;
`grid_container_ext` provides hooks for optional modules (multigrid etc.) to store their data in `cg`.
Fine/coarse flux helpers: `fc_fluxes`, `flx_arr`, `flx_cell`, `grid_container_fcflx`.

## Lists and levels
Chain: `cg_list` → `cg_list_misc` → `cg_list_dataop` (field arithmetic, norms) → `cg_list_bnd` (internal/external boundaries, MPI
exchange) → `cg_list_balance` → `cg_list_neighbors` → `cg_level` → `cg_level_connected` (inter-level prolong/restrict, vertical
boundary prep; the biggest file) → `cg_level_base`, `cg_level_coarsest`, `cg_level_finest` (specialised levels).
`cg_leaves` (all finest blocks – the standard iteration target `leaves`), `cg_list_global` (`all_cg`), `list_of_cg_lists`
(`all_lists`), `grid.F90` (init/cleanup of the above), `find_level`, `gcpa`, `tag_pool`, `merge_segments`.

## Load balancing
`load_balance`, `load_balance_helpers` (`lb_helpers`), `rebalance`, `cg_cost*` (cost accounting), `cg_list_balance`, `sort_*`
(sortable lists of cg, pieces, segments – built on `base/sortable_list`).

## AMR
`refinement.F90` (parameters, namelist `AMR`), `refinement_flag` (flags), `refinement_update` (refine/derefine driver, called from the main loop),
`cg_expand_base`. Unified Refinement Criteria are in `refinement/`: `URC` base type, `URC_list`, `URC_var` (field-based filters),
`URC_Jeans`, `URC_nbody`, `URC_user`, `URC_filter`, `URCgeom*` (box, point, z-cylinder). New criterion = new `URC_*` type + registration in `URC_list`.

## Gotchas
- Iteration idiom: `cgl => leaves%first; do while (associated(cgl)); cg => cgl%cg; ...; cgl => cgl%nxt; enddo`.
- Do not reorder `use` statements into cycles: several modules (`level_essentials`, `cg_list_*` chain) exist solely to break circular dependencies.
- Refinement tests: `problems/sedov/problem.par.AMR`, gold test `MHDsedovAMR` (`make gold_MHDsedovAMR`) – run them after touching `cg_level_connected`, `refinement_update`, `cg_list_bnd`.
