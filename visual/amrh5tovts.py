#!/usr/bin/env python3
import sys
import h5py
import numpy as np

# --- VTK Imports ---
try:
    from vtk import (
        vtkOverlappingAMR,
        vtkAMRBox,
        vtkUniformGrid,
        vtkCellData,
        vtkFloatArray,
        vtkXMLUniformGridAMRWriter,   # NEW WRITER
        VTK_XYZ_GRID, VTK_XY_PLANE, VTK_XZ_PLANE, VTK_YZ_PLANE
    )
    from vtk.util import numpy_support
except ImportError:
    print("Error: The 'vtk' Python package is required.")
    print("Please install it, for example using: pip install vtk")
    sys.exit(1)


def get_field_names(h5_file):
    print("Scanning for available data fields...")
    first_block_name = list(h5_file["data"].keys())[0]
    field_names = list(h5_file["data"][first_block_name].keys())
    print(f"Found fields: {field_names}")
    return field_names


def create_vtk_amr_dataset(h5_file, field_names):
    # --- 1) Read block metadata ---
    print("Reading block metadata...")
    block_meta = []
    all_levels = set()
    for name in h5_file["data"].keys():
        g = h5_file["data"][name]
        level = int(np.array(g.attrs["level"]).ravel()[0])
        all_levels.add(level)
        meta = {
            "name": name,
            "level": level,
            "offset": g.attrs["off"].astype(int),
            "dims":   g.attrs["n_b"].astype(int),
            "spacing": g.attrs["dl"],
            "origin":  g.attrs["left_edge"],
        }
        block_meta.append(meta)

    # Determine if dataset is 2-D (one axis==1 for all blocks)
    collapsed_axis = None
    for ax in range(3):
        if all(m["dims"][ax] == 1 for m in block_meta):
            collapsed_axis = ax
            break

    # --- 2) Init AMR ---
    print("Initializing VTK AMR structure...")
    num_levels = len(all_levels)
    sim_params = h5_file["simulation_parameters"].attrs
    refine_by = int(np.array(sim_params["refine_by"]).ravel()[0])

    coarsest_meta = next((m for m in block_meta if m["level"] == 0), None)
    if not coarsest_meta:
        raise ValueError("No level-0 blocks found.")

    origin = sim_params["domain_left_edge"]
    coarsest_spacing = coarsest_meta["spacing"]

    blocks_per_level = [sum(1 for m in block_meta if m["level"] == i) for i in range(num_levels)]

    amr_dataset = vtkOverlappingAMR()
    amr_dataset.Initialize(num_levels, blocks_per_level)

    # Grid description: true 2-D when possible
    if collapsed_axis is None:
        amr_dataset.SetGridDescription(VTK_XYZ_GRID)
    else:
        plane = {0: VTK_YZ_PLANE, 1: VTK_XZ_PLANE, 2: VTK_XY_PLANE}[collapsed_axis]
        amr_dataset.SetGridDescription(plane)

    amr_dataset.SetOrigin(origin)
    for i in range(num_levels):
        level_spacing = coarsest_spacing / (refine_by**i)
        amr_dataset.SetSpacing(i, level_spacing)

    # --- 3) Fill AMR ---
    print("Defining AMR boxes and loading data for each block...")
    block_meta.sort(key=lambda m: (m["level"], m["name"]))
    block_idx_per_level = [0] * num_levels

    for meta in block_meta:
        level = meta["level"]
        offset = meta["offset"]
        dims = meta["dims"].copy()  # cells (nx, ny, nz)

        # AMR box in cell indices
        lo_corner = offset
        hi_corner = offset + dims - 1
        box = vtkAMRBox(lo_corner, hi_corner)

        block_idx = block_idx_per_level[level]
        amr_dataset.SetAMRBox(level, block_idx, box)

        grid = vtkUniformGrid()
        grid.SetOrigin(meta["origin"])
        grid.SetSpacing(meta["spacing"])

        # POINT dims need special handling for 2-D:
        pdims = dims + 1  # point counts
        if collapsed_axis is not None:
            pdims[collapsed_axis] = 1  # true 2-D plane (no 1-cell thickness)

        grid.SetDimensions(int(pdims[0]), int(pdims[1]), int(pdims[2]))

        # Cell-centered arrays
        cell_data = grid.GetCellData()
        for field in field_names:
            data_zyx = h5_file["data"][meta["name"]][field][:]
            flat_array = np.asarray(data_zyx, dtype=np.float32).ravel(order="C")
            vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
            vtk_array.SetName(field)
            cell_data.AddArray(vtk_array)

        amr_dataset.SetDataSet(level, block_idx, grid)
        block_idx_per_level[level] += 1

    print("Finished loading all blocks.")
    return amr_dataset


def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_piernik_amr_to_vtk.py <input_file.h5> <output_file.vth|.vthb>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    print(f"Starting conversion of '{input_filename}' to '{output_filename}'")

    with h5py.File(input_filename, "r") as f:
        field_names = get_field_names(f)
        vtk_amr_dataset = create_vtk_amr_dataset(f, field_names)

    print(f"Writing VTK file: {output_filename}")
    writer = vtkXMLUniformGridAMRWriter()
    writer.SetFileName(output_filename)
    writer.SetInputData(vtk_amr_dataset)

# IMPORTANT: write the meta/root file (.vth or .vthb)
    if hasattr(writer, "SetWriteMetaFile"):
       writer.SetWriteMetaFile(1)   # <-- this is the key change

# Optional: keep piece files in the same directory (no extra subfolder)
    if hasattr(writer, "SetUseSubdirectory"):
      writer.SetUseSubdirectory(0)

# If you pass a .vthb filename, try to pack data into the root (single-file if supported)
    ext = output_filename.lower()
    if ext.endswith(".vthb"):
        if hasattr(writer, "SetDataModeToAppended"): writer.SetDataModeToAppended()
        if hasattr(writer, "EncodeAppendedDataOff"): writer.EncodeAppendedDataOff()
        if hasattr(writer, "SetCompressorTypeToNone"): writer.SetCompressorTypeToNone()

    ok = writer.Write()
    if ok == 0:
       raise RuntimeError("VTK writer failed (no file written).")

    print("\nConversion complete!")


if __name__ == "__main__":
    main()

