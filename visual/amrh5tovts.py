#!/usr/bin/env python3
"""
Converts Piernik-style AMR (Adaptive Mesh Refinement) HDF5 simulation data
into a single VTK Overlapping AMR file (.vthb or .vth) suitable for
visualization in tools like ParaView.

The script correctly handles both 2D (as a thin plane) and 3D datasets,
calculates the necessary AMR metadata, and uses the modern VTK writer to
produce a single, portable output file where possible.
"""

import sys
import h5py
import numpy as np
import multiprocessing as mp
import argparse
import os
from functools import partial
# --- VTK Imports ---
try:
    from vtk import (
        vtkOverlappingAMR,
        vtkAMRBox,
        vtkUniformGrid,
        vtkXMLUniformGridAMRWriter,
        VTK_XYZ_GRID, VTK_XY_PLANE, VTK_XZ_PLANE, VTK_YZ_PLANE
    )
    from vtk.util import numpy_support
except ImportError:
    print("Error: The 'vtk' Python package is required.", file=sys.stderr)
    print("Please install it, for example using: pip install vtk", file=sys.stderr)
    sys.exit(1)


def get_field_names(h5_file):
    """
    Scans the first data block in the HDF5 file to find all field names.

    Args:
        h5_file (h5py.File): An open HDF5 file object.

    Returns:
        list[str]: A list of the dataset names (fields).
    """
    print("Scanning for available data fields...")
    try:
        first_block_name = next(iter(h5_file["data"]))
        field_names = list(h5_file["data"][first_block_name].keys())
        print(f"Found fields: {field_names}")
        return field_names
    except (StopIteration, KeyError):
        raise ValueError("HDF5 file does not contain a '/data' group or it is empty.")


def create_vtk_amr_dataset(h5_file, field_names):
    """
    Builds the vtkOverlappingAMR dataset from the HDF5 source file.

    Args:
        h5_file (h5py.File): The open HDF5 file.
        field_names (list[str]): The names of the fields to process.

    Returns:
        vtkOverlappingAMR: The fully constructed VTK AMR dataset.
    """
    # --- 1. Gather Block Metadata ---
    print("Reading block metadata...")
    blocks = []
    levels_set = set()
    for name, g in h5_file["data"].items():
        lvl = int(np.array(g.attrs["level"]).ravel()[0])
        levels_set.add(lvl)
        blocks.append(dict(
            name=name,
            level=lvl,
            dims=g.attrs["n_b"].astype(int),
            spacing=np.array(g.attrs["dl"], float),
            origin=np.array(g.attrs["left_edge"], float)
        ))

    levels = sorted(levels_set)
    level_to_index = {L: i for i, L in enumerate(levels)}
    num_levels = len(levels)

    # --- 2. Determine Dimensionality and Per-Level Spacing ---
    collapsed_axis = None
    for ax in range(3):
        if all(b["dims"][ax] == 1 for b in blocks):
            collapsed_axis = ax
            break

    level_spacing = {
        L: next(b["spacing"] for b in blocks if b["level"] == L)
        for L in levels
    }

    # --- 3. Initialize AMR Container ---
    print("Initializing VTK AMR container...")
    amr = vtkOverlappingAMR()
    blocks_per_level = [sum(1 for b in blocks if b["level"] == L) for L in levels]
    amr.Initialize(num_levels, blocks_per_level)

    if collapsed_axis is None:
        amr.SetGridDescription(VTK_XYZ_GRID)
    else:
        plane_map = {0: VTK_YZ_PLANE, 1: VTK_XZ_PLANE, 2: VTK_XY_PLANE}
        amr.SetGridDescription(plane_map[collapsed_axis])

    try:
        domain_left = np.array(h5_file["simulation_parameters"].attrs["domain_left_edge"], float)
    except KeyError:
        domain_left = np.min([b["origin"] for b in blocks], axis=0)
    amr.SetOrigin(domain_left)

    for i, L in enumerate(levels):
        amr.SetSpacing(i, level_spacing[L])
        if i > 0:
            prev_level = levels[i - 1]
            ratio = level_spacing[prev_level][0] / level_spacing[L][0]
            amr.SetRefinementRatio(i, int(round(ratio)))

    # --- 4. Populate AMR with Boxes and Grid Data ---
    print("Populating AMR with block data...")
    blocks.sort(key=lambda b: (b["level"], b["name"]))
    next_idx = [0] * num_levels

    for b in blocks:
        level_idx = level_to_index[b["level"]]
        block_idx = next_idx[level_idx]
        dims = b["dims"]

        # AMRBox is defined in the index space of its own level.
        lo = np.rint((b["origin"] - domain_left) / level_spacing[b["level"]]).astype(int)
        hi = lo + dims - 1
        box = vtkAMRBox(lo, hi)
        amr.SetAMRBox(level_idx, block_idx, box)

        # Create the vtkUniformGrid for this block.
        ug = vtkUniformGrid()
        ug.SetOrigin(b["origin"])
        ug.SetSpacing(b["spacing"])

        point_dims = dims + 1
        if collapsed_axis is not None:
            point_dims[collapsed_axis] = 1
        ug.SetDimensions(int(point_dims[0]), int(point_dims[1]), int(point_dims[2]))

        # Attach cell data arrays.
        cell_data = ug.GetCellData()
        for field in field_names:
            arr_data = h5_file["data"][b["name"]][field][:]
            flat_arr = arr_data.astype(np.float32).ravel(order="C")
            vtk_arr = numpy_support.numpy_to_vtk(flat_arr, deep=True)
            vtk_arr.SetName(field)
            cell_data.AddArray(vtk_arr)

        amr.SetDataSet(level_idx, block_idx, ug)
        next_idx[level_idx] += 1

    for i, L in enumerate(levels):
        print(f"[build] Level {L} (Index {i}): {blocks_per_level[i]} blocks, spacing={level_spacing[L]}")

    return amr



def convert_file(input_path, output_dir):
    base = os.path.basename(input_path).replace(".h5", ".vthb")
    output_path = os.path.join(output_dir, base)

    try:
        with h5py.File(input_path, "r") as f:
            field_names = get_field_names(f)
            vtk_amr_dataset = create_vtk_amr_dataset(f, field_names)

        print(f"[write] Writing VTK file: {output_path}")
        writer = vtkXMLUniformGridAMRWriter()
        writer.SetFileName(output_path)
        writer.SetInputData(vtk_amr_dataset)

        if output_path.lower().endswith(".vthb"):
            if hasattr(writer, "SetDataModeToAppended"):
                writer.SetDataModeToAppended()
            if hasattr(writer, "EncodeAppendedDataOff"):
                writer.EncodeAppendedDataOff()

        if hasattr(writer, "SetUseSubdirectory"):
            writer.SetUseSubdirectory(0)

        if writer.Write() == 0:
            raise RuntimeError("VTK writer failed.")

        print(f"[done] {input_path} → {output_path}")
    except Exception as e:
        print(f"[error] Failed to convert {input_path}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Convert Piernik AMR HDF5 to .vthb (VTK AMR)")
    parser.add_argument("input", nargs="?", help="Input .h5 file (if omitted, all *.h5 files in directory)")
    parser.add_argument("output", nargs="?", help="Output .vthb file (only with single input)")
    parser.add_argument("-j", "--jobs", type=int, default=mp.cpu_count(),
                        help="Number of parallel jobs [default: all cores]")
    args = parser.parse_args()

    output_dir = "vthb_out"
    os.makedirs(output_dir, exist_ok=True)

    if args.input is None:
        # Batch mode: all .h5 files
        files = sorted(f for f in os.listdir(".") if f.endswith(".h5"))
        print(f"[info] Converting {len(files)} .h5 files using {args.jobs} processes...")

        with mp.Pool(processes=args.jobs) as pool:
            pool.map(partial(convert_file, output_dir=output_dir), files)
    else:
        if args.output is None:
            base = os.path.splitext(os.path.basename(args.input))[0]
            output = os.path.join(output_dir, base + ".vthb")
        else:
            output = os.path.join(output_dir, os.path.basename(args.output))
        convert_file(args.input, output_dir)

if __name__ == "__main__":
    main()



