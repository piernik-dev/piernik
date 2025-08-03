#!/usr/bin/env python3
#
#  Python script to convert Piernik-style multi-block HDF5 files with AMR
#  to a single VTK Overlapping AMR (.vth) file for ParaView.
#
#  This script analyzes the HDF5 file structure, describes the multi-level
#  block hierarchy to VTK, and saves all scalar and vector fields as
#  CELL-CENTERED data.
#
#  Author: Gemini
#  Date: August 2, 2025
#
#  Usage:
#     python convert_piernik_amr_to_vtk.py <input_file.h5> <output_file.vth>
#
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
        vtkXMLHierarchicalBoxDataWriter,
        VTK_XYZ_GRID
    )
    from vtk.util import numpy_support
except ImportError:
    print("Error: The 'vtk' Python package is required.")
    print("Please install it, for example using: pip install vtk")
    sys.exit(1)


def get_field_names(h5_file):
    """
    Scans the first data block in the HDF5 file to find all available
    dataset names.
    """
    print("Scanning for available data fields...")
    field_names = []
    first_block_name = list(h5_file["data"].keys())[0]
    for name in h5_file["data"][first_block_name]:
        field_names.append(name)
    print(f"Found fields: {field_names}")
    return field_names


def create_vtk_amr_dataset(h5_file, field_names):
    """
    Creates a vtkOverlappingAMR dataset from the HDF5 file.

    Args:
        h5_file (h5py.File): An open H5py file object.
        field_names (list): A list of strings with the names of the data fields.

    Returns:
        vtkOverlappingAMR: The fully populated VTK AMR dataset.
    """
    # --- 1. Read all block metadata first ---
    print("Reading block metadata...")
    block_meta = []
    all_levels = set()
    for name in h5_file["data"].keys():
        block_group = h5_file["data"][name]
        level = int(block_group.attrs["level"][0])
        all_levels.add(level)
        meta = {
            "name": name,
            "level": level,
            "offset": block_group.attrs["off"].astype(int),
            "dims": block_group.attrs["n_b"].astype(int),
            "spacing": block_group.attrs["dl"],
            "origin": block_group.attrs["left_edge"],
        }
        block_meta.append(meta)

    # --- 2. Initialize the VTK AMR dataset ---
    print("Initializing VTK AMR structure...")
    num_levels = len(all_levels)
    sim_params = h5_file["simulation_parameters"].attrs

    refine_by = int(sim_params["refine_by"][0])

    coarsest_meta = next((m for m in block_meta if m["level"] == 0), None)
    if not coarsest_meta:
        raise ValueError("Could not find any blocks at the coarsest level (level 0).")

    origin = sim_params["domain_left_edge"]
    coarsest_spacing = coarsest_meta["spacing"]

    # --- FIX: Pre-calculate the number of blocks per level ---
    blocks_per_level = [sum(1 for m in block_meta if m["level"] == i) for i in range(num_levels)]

    amr_dataset = vtkOverlappingAMR()

    # --- FIX: Initialize the AMR object with the number of levels AND the block counts ---
    amr_dataset.Initialize(num_levels, blocks_per_level)

    amr_dataset.SetGridDescription(VTK_XYZ_GRID)
    amr_dataset.SetOrigin(origin)

    for i in range(num_levels):
        level_spacing = coarsest_spacing / (refine_by**i)
        amr_dataset.SetSpacing(i, level_spacing)

    # --- 3. Define the AMR block structure and load data ---
    print("Defining AMR boxes and loading data for each block...")

    block_meta.sort(key=lambda m: (m["level"], m["name"]))

    block_idx_per_level = [0] * num_levels
    for meta in block_meta:
        level = meta["level"]
        offset = meta["offset"]
        dims = meta["dims"]

        lo_corner = offset
        hi_corner = offset + dims - 1

        box = vtkAMRBox(lo_corner, hi_corner)

        block_idx = block_idx_per_level[level]
        amr_dataset.SetAMRBox(level, block_idx, box)

        grid = vtkUniformGrid()
        grid.SetOrigin(meta["origin"])
        grid.SetSpacing(meta["spacing"])
        grid.SetDimensions(dims[0] + 1, dims[1] + 1, dims[2] + 1)

        cell_data = grid.GetCellData()
        for field in field_names:
            data_zyx = h5_file["data"][meta["name"]][field][:]
            flat_array = data_zyx.flatten(order="C")

            vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
            vtk_array.SetName(field)
            cell_data.AddArray(vtk_array)

        amr_dataset.SetDataSet(level, block_idx, grid)

        block_idx_per_level[level] += 1

    print("Finished loading all blocks.")
    return amr_dataset


def main():
    """
    Main function to drive the conversion process.
    """
    if len(sys.argv) != 3:
        print("Usage: python convert_piernik_amr_to_vtk.py <input_file.h5> <output_file.vth>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    print(f"Starting conversion of '{input_filename}' to '{output_filename}'")

    with h5py.File(input_filename, "r") as f:
        field_names = get_field_names(f)
        vtk_amr_dataset = create_vtk_amr_dataset(f, field_names)

    print(f"Writing VTK file: {output_filename}")
    writer = vtkXMLHierarchicalBoxDataWriter()
    writer.SetFileName(output_filename)
    writer.SetInputData(vtk_amr_dataset)
    writer.Write()

    print("\nConversion complete!")


if __name__ == "__main__":
    main()
