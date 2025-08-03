#!/usr/bin/env python3
#
#  Python script to convert Piernik-style multi-block HDF5 files
#  to a single VTK Structured Grid (.vts) file for ParaView.
#
#  This script analyzes the HDF5 file structure, stitches together
#  all the data blocks into a single contiguous grid, and saves all
#  scalar and vector fields as CELL-CENTERED data.
#
#  Author: Gemini
#  Date: August 1, 2025
#
#  Usage:
#     python convert_piernik_to_vtk_cell_centered.py <input_file.h5> <output_file.vts>
#

import sys
import h5py
import numpy as np

# --- VTK Imports ---
# We wrap this in a try/except block to provide a helpful error message
# if the user doesn't have VTK installed.
try:
    from vtk import (
        vtkStructuredGrid,
        vtkPoints,
        vtkFloatArray,
        vtkXMLStructuredGridWriter,
    )
    from vtk.util import numpy_support
except ImportError:
    print("Error: The 'vtk' Python package is required.")
    print("Please install it, for example using: pip install vtk")
    sys.exit(1)


def get_field_names(h5_file):
    """
    Scans the first data block in the HDF5 file to find all available
    dataset names (e.g., 'density', 'pressure').

    Args:
        h5_file (h5py.File): An open H5py file object.

    Returns:
        list: A list of strings with the names of the data fields.
    """
    print("Scanning for available data fields...")
    field_names = []
    # Get the name of the first grid block (e.g., 'grid_0000000000')
    first_block_name = list(h5_file["data"].keys())[0]
    for name in h5_file["data"][first_block_name]:
        field_names.append(name)
    print(f"Found fields: {field_names}")
    return field_names


def load_and_stitch_data(fname, field_names):
    """
    Loads and stitches a 3D field from multiple blocks in an HDF5 file.
    This function is based on the user-provided reference code and the
    HDF5 inspection report.

    Args:
        fname (str): Path to the HDF5 file.
        field_names (list): A list of all field names to extract.

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary of the stitched 3D data arrays for each field.
            - np.ndarray: The global dimensions of the stitched cell grid (Nx, Ny, Nz).
            - np.ndarray: The physical coordinates of the grid origin (x0, y0, z0).
            - np.ndarray: The cell spacing in each direction (dx, dy, dz).
    """
    stitched_data = {}

    with h5py.File(fname, "r") as f:
        # --- 1. Read metadata from the HDF5 file ---
        block_names = list(f["data"].keys())
        all_offsets = np.array([f["data"][name].attrs["off"] for name in block_names], dtype=int)
        all_dims = np.array([f["data"][name].attrs["n_b"] for name in block_names], dtype=int)

        # Determine the full size of the stitched grid of cells
        global_cell_dims = np.max(all_offsets + all_dims, axis=0)
        print(f"Global cell grid dimensions determined to be: {global_cell_dims}")

        # Get global grid physical properties from the root attributes
        origin = f["simulation_parameters"].attrs["domain_left_edge"]
        domain_size = f["simulation_parameters"].attrs["domain_right_edge"] - origin

        # For cell-centered data, spacing is domain size / number of cells.
        spacing = domain_size / np.maximum(1, global_cell_dims)

        print(f"Global grid origin: {origin}")
        print(f"Cell spacing: {spacing}")

        # --- 2. Load and stitch each data field ---
        for field in field_names:
            print(f"  Processing field: {field}")
            # Create an empty array to hold the full stitched data
            # The data on disk is (z, y, x), so we create the numpy array in that order.
            stitched_field = np.empty((global_cell_dims[2], global_cell_dims[1], global_cell_dims[0]), dtype=np.float32)

            # Iterate over each block and place its data into the large array
            for i, name in enumerate(block_names):
                block_data = f["data"][name][field][:]  # This is already (z, y, x)
                offset = all_offsets[i]
                dims = all_dims[i]

                # Define the slice where this block goes in the big array
                slc = np.s_[
                    offset[2]: offset[2] + dims[2],
                    offset[1]: offset[1] + dims[1],
                    offset[0]: offset[0] + dims[0],
                ]
                stitched_field[slc] = block_data

            stitched_data[field] = stitched_field

    return stitched_data, global_cell_dims.astype(int), origin, spacing


def create_vtk_grid(stitched_data, cell_dims, origin, spacing):
    """
    Creates a vtkStructuredGrid object from the stitched data and metadata,
    storing the data as cell-centered.

    Args:
        stitched_data (dict): Dictionary of stitched NumPy arrays for each field.
        cell_dims (np.ndarray): The global dimensions of the cell grid (Nx, Ny, Nz).
        origin (np.ndarray): The physical coordinates of the grid origin.
        spacing (np.ndarray): The cell spacing in each direction.

    Returns:
        vtkStructuredGrid: The fully populated VTK grid object.
    """
    # --- 1. Create the vtkStructuredGrid object and set its dimensions ---
    grid = vtkStructuredGrid()

    # --- MODIFIED LOGIC ---
    # Determine the number of points (nodes) needed to define the grid.
    # If a dimension has only 1 cell (e.g., a 2D plane), it's represented by 1 layer of points.
    # If a dimension has N>1 cells, it needs N+1 layers of points.
    # An edge case is a single cell [1,1,1], which requires [2,2,2] points.
    if np.all(cell_dims == 1):
        # Handle the edge case of a single 3D cell
        point_dims = cell_dims + 1
    else:
        # For 1D, 2D, or 3D grids, calculate point dimensions.
        # A dimension with >1 cell needs N+1 points.
        # A dimension with 1 cell is collapsed into a single layer of points (dim=1).
        point_dims = np.copy(cell_dims)
        point_dims[cell_dims > 1] += 1

    # The dimensions of a structured grid in VTK are the number of points.
    grid.SetDimensions(point_dims[0], point_dims[1], point_dims[2])

    # --- 2. Create the physical coordinates of the grid points (nodes) ---
    points = vtkPoints()
    points.SetNumberOfPoints(np.prod(point_dims))

    print("Generating grid corner points (nodes)...")
    # Generate coordinates for the points. If a point dimension is 1, np.arange(1)
    # correctly creates a single coordinate for that axis.
    x_coords = origin[0] + spacing[0] * np.arange(point_dims[0])
    y_coords = origin[1] + spacing[1] * np.arange(point_dims[1])
    z_coords = origin[2] + spacing[2] * np.arange(point_dims[2])

    # Create a meshgrid and reshape it to a list of (x,y,z) points.
    zz, yy, xx = np.meshgrid(z_coords, y_coords, x_coords, indexing='ij')
    point_coords = np.vstack([xx.ravel('C'), yy.ravel('C'), zz.ravel('C')]).T

    # Convert the numpy array of points to a VTK points array
    vtk_points_array = numpy_support.numpy_to_vtk(point_coords.astype(np.float32), deep=True)
    points.SetData(vtk_points_array)
    grid.SetPoints(points)

    # --- 3. Add data arrays (scalars and vectors) to the grid's CELLS ---
    print("Adding data arrays to VTK grid cells...")
    cell_data = grid.GetCellData()

    for field_name, data_array in stitched_data.items():
        # The data array is (z, y, x), matching the cell layout.
        # Flatten it to a 1D array with "C" order (last index varies fastest).
        flat_array = data_array.flatten(order="C")

        # Convert the NumPy array to a VTK array
        vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
        vtk_array.SetName(field_name)
        cell_data.AddArray(vtk_array)
        print(f"  Added cell-centered scalar field: {field_name}")

    return grid


def main():
    """
    Main function to drive the conversion process.
    """
    if len(sys.argv) != 3:
        print("Usage: python convert_piernik_to_vtk_cell_centered.py <input_file.h5> <output_file.vts>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    print(f"Starting conversion of '{input_filename}' to '{output_filename}'")

    # --- Step 1: Discover what data fields are in the file ---
    with h5py.File(input_filename, "r") as f:
        field_names = get_field_names(f)

    # --- Step 2: Load all fields and stitch them into global arrays ---
    stitched_data, cell_dims, origin, spacing = load_and_stitch_data(
        input_filename, field_names
    )

    # --- Step 3: Create a VTK grid object from the stitched data ---
    vtk_grid = create_vtk_grid(stitched_data, cell_dims, origin, spacing)

    # --- Step 4: Write the VTK grid to a .vts file ---
    print(f"Writing VTK file: {output_filename}")
    writer = vtkXMLStructuredGridWriter()
    writer.SetFileName(output_filename)
    writer.SetInputData(vtk_grid)
    writer.Write()

    print("\nConversion complete!")


if __name__ == "__main__":
    main()
