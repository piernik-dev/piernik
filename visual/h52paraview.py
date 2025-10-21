#!/usr/bin/env python3


"""
Converts PIERNIK output .h5 file to Paraview friendly .vthb file.
Handles both AMR and non-AMR dataset.
Arguments : -i [input_file.h5/no arguments] -o [output_file.vthd/no arguments] -j [nprocs] -v [verbosity level:0/1]
If no arguments are given for input and output then all files inside the parent directory will be converted
Note that all converted files are stored in
"""

import sys
import h5py
import numpy as np
import multiprocessing as mp
import argparse
import os
from functools import partial

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
    print("To install use: pip install vtk", file=sys.stderr)
    sys.exit(1)


def get_field_names(h5_file):
    print("[INFO] Scanning for available data fields...")
    try:
        first_block_name = next(iter(h5_file["data"]))
        field_names = list(h5_file["data"][first_block_name].keys())
        print(f"[INFO : ] Found fields: {field_names}")
        return field_names
    except (StopIteration, KeyError):
        raise ValueError("HDF5 file does not contain a '/data' group or it is empty.")


def create_vtk_amr_dataset(h5_file, field_names, v):
    if v == 1:
        print("[INFO] Reading block metadata...")
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
    collapsed_axis = None
    for ax in range(3):
        if all(b["dims"][ax] == 1 for b in blocks):
            collapsed_axis = ax
            break
    level_spacing = {
        L: next(b["spacing"] for b in blocks if b["level"] == L)
        for L in levels
    }
    if v == 1:
        print("[INFO] Initializing VTK AMR container...")
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
    if v == 1:
        print("[INFO] Populating AMR with block data...")
    blocks.sort(key=lambda b: (b["level"], b["name"]))
    next_idx = [0] * num_levels
    for b in blocks:
        level_idx = level_to_index[b["level"]]
        block_idx = next_idx[level_idx]
        dims = b["dims"]
        lo = np.rint((b["origin"] - domain_left) / level_spacing[b["level"]]).astype(int)
        hi = lo + dims - 1
        box = vtkAMRBox(lo, hi)
        amr.SetAMRBox(level_idx, block_idx, box)
        ug = vtkUniformGrid()
        ug.SetOrigin(b["origin"])
        ug.SetSpacing(b["spacing"])
        point_dims = dims + 1
        if collapsed_axis is not None:
            point_dims[collapsed_axis] = 1
        ug.SetDimensions(int(point_dims[0]), int(point_dims[1]), int(point_dims[2]))
        cell_data = ug.GetCellData()
        for field in field_names:
            arr_data = h5_file["data"][b["name"]][field][:]
            flat_arr = arr_data.astype(np.float32).ravel(order="C")
            vtk_arr = numpy_support.numpy_to_vtk(flat_arr, deep=True)
            vtk_arr.SetName(field)
            cell_data.AddArray(vtk_arr)
        amr.SetDataSet(level_idx, block_idx, ug)
        next_idx[level_idx] += 1
    if v == 1:
        for i, L in enumerate(levels):
            print(f"[BUILD] Level {L} (Index {i}): {blocks_per_level[i]} blocks, spacing={level_spacing[L]}")
    return amr


def convert_file(input_path, output_dir, v, output_path=None):
    if output_path is None:
        base = os.path.basename(input_path).replace(".h5", ".vthb")
        output_path = os.path.join(output_dir, base)
    nlevels_detected = 1
    try:
        with h5py.File(input_path, "r") as f:
            field_names = get_field_names(f)
            vtk_amr_dataset = create_vtk_amr_dataset(f, field_names, v=v)
            try:
                nlevels_detected = int(vtk_amr_dataset.GetNumberOfLevels())
            except Exception:
                nlevels_detected = 1
        print(f"[WRITE] Writing VTK file: {output_path}")
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
        print(f"[DONE] {input_path} → {output_path}")
    except Exception as e:
        print(f"[ERROR] Failed to convert {input_path}: {e}")
    return nlevels_detected


def main():
    parser = argparse.ArgumentParser(description="Convert Piernik AMR HDF5 to .vthb (VTK AMR) compatible with Paraview")
    parser.add_argument("input", nargs="?", help="Input .h5 file (if omitted, all *.h5 files in directory)")
    parser.add_argument("output", nargs="?", help="Output .vthb file (only with single input)")
    parser.add_argument("-i", "--input", dest="input_opt", help="Same as positional input")
    parser.add_argument("-o", "--output", dest="output_opt", help="Same as positional output")
    parser.add_argument("-j", "--jobs", type=int, default=mp.cpu_count(),
                        help="Number of parallel jobs [default: all cores]")
    parser.add_argument(
        "-v", "--verbosity",
        metavar="LEVEL",
        type=int,
        choices=[0, 1],
        default=0,
        help="Verbosity level 0 (brief) or 1 (detailed) [default: 0]"
    )
    args = parser.parse_args()
    infile = args.input or args.input_opt
    outfile = args.output or args.output_opt
    output_dir = "vthb_out"
    os.makedirs(output_dir, exist_ok=True)
    if infile is None:
        files = sorted(f for f in os.listdir(".") if f.endswith(".h5"))
        if not files:
            print("[ERROR] No input provided and no *.h5 files found in the current directory.", file=sys.stderr)
            sys.exit(2)
        if outfile is not None:
            print("[ERROR] --output/OUTPUT is only valid with a single input file.", file=sys.stderr)
            sys.exit(2)
        print(f"[INFO] Converting {len(files)} .h5 files using {args.jobs} processes...")
        if args.jobs == 1:
            levels = []
            for f in files:
                levels.append(convert_file(f, output_dir, v=args.verbosity))
        else:
            with mp.Pool(processes=args.jobs) as pool:
                levels = pool.starmap(convert_file, [(f, output_dir, args.verbosity) for f in files])
        if any(n > 1 for n in levels):
            print("[PRO TIP] Detected AMR with multiple levels. In ParaView, increase "
                  "the ‘Default Number of Levels’ or set it 0 in the Information tab to see all levels.")
    else:
        if outfile is not None:
            out_name = os.path.basename(outfile)
            explicit_out = os.path.join(output_dir, out_name)
        else:
            explicit_out = None
        nlevels = convert_file(infile, output_dir, v=args.verbosity, output_path=explicit_out)
        if nlevels > 1:
            print("[PRO TIP] Detected AMR with multiple levels. In ParaView, increase "
                  "the ‘Default Number of Levels’ or set it 0 in the Information tab to see all levels.")


if __name__ == "__main__":
    main()
