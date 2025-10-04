#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit


def load_piernik_data(fname, field_name="density"):
    """
    Reads and stitches together Piernik's block-based HDF5 files,
    correctly reconstructing physical coordinates from per-block attributes.
    """
    print(f"INFO: Loading '{field_name}' from '{fname}'...")
    with h5py.File(fname, "r") as f:
        # --- START OF FIX ---
        # Read the index/dimension datasets. They are stored with (x, y, z) columns.
        starts_xyz = f["grid_left_index"][:]
        dims_xyz = f["grid_dimensions"][:]

        permutation_zyx = [2, 1, 0]
        starts_zyx = starts_xyz[:, permutation_zyx]
        dims_zyx = dims_xyz[:, permutation_zyx]
        # --- END OF FIX ---

        global_shape_zyx = np.max(starts_zyx + dims_zyx, axis=0)

        data = np.empty(global_shape_zyx, dtype=np.float64)
        edges_zyx = [np.empty(g + 1, dtype=np.float64) for g in global_shape_zyx]

        for i, block_name in enumerate(f["data"]):
            block_group = f["data"][block_name]
            block_data = block_group[field_name][:]

            slc_zyx = tuple(slice(starts_zyx[i, d], starts_zyx[i, d] + dims_zyx[i, d]) for d in range(3))

            data[slc_zyx] = block_data

            if 'left_edge' in block_group.attrs and 'dl' in block_group.attrs:
                left_edge_xyz = block_group.attrs['left_edge']
                dcell_xyz = block_group.attrs['dl']

                n_z, i0_z = dims_zyx[i, 0], starts_zyx[i, 0]
                loc_edges_z = left_edge_xyz[2] + dcell_xyz[2] * np.arange(n_z + 1)
                edges_zyx[0][i0_z:i0_z + n_z + 1] = loc_edges_z

                n_y, i0_y = dims_zyx[i, 1], starts_zyx[i, 1]
                loc_edges_y = left_edge_xyz[1] + dcell_xyz[1] * np.arange(n_y + 1)
                edges_zyx[1][i0_y:i0_y + n_y + 1] = loc_edges_y

                n_x, i0_x = dims_zyx[i, 2], starts_zyx[i, 2]
                loc_edges_x = left_edge_xyz[0] + dcell_xyz[0] * np.arange(n_x + 1)
                edges_zyx[2][i0_x:i0_x + n_x + 1] = loc_edges_x

        zc = 0.5 * (edges_zyx[0][:-1] + edges_zyx[0][1:])
        yc = 0.5 * (edges_zyx[1][:-1] + edges_zyx[1][1:])
        xc = 0.5 * (edges_zyx[2][:-1] + edges_zyx[2][1:])

    return data, zc, yc, xc


# --- 2. Function Definition for the Sine Wave Model ---
def sine_wave_model(x, amplitude, wavenumber, phase, offset):
    """Model function for a sine wave: f(x) = A*sin(k*x + phi) + D"""
    return amplitude * np.sin(wavenumber * x + phase) + offset


file_initial = 'cpaw_tst_0000.h5'
file_final = 'cpaw_tst_0001.h5'                # Make sure you choose the correct file coresponding to t = T
variables_to_plot = ['mag_field_y', 'mag_field_z']

# Theoretical wavenumbers from your parameter file (for initial guesses)
kx_theory = 2 * np.pi
ky_theory = 4 * np.pi

# --- Create the 2x2 figure BEFORE the loop ---
fig, axes = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)
fig.suptitle(r'1D Slice of Magnetic Field Components Amp * $\sin(kx + \delta) + D$', fontsize=12)

for i, var_name in enumerate(variables_to_plot):
    try:
        data_init, zc, yc, xc = load_piernik_data(file_initial, field_name=var_name)
        data_final, _, _, _ = load_piernik_data(file_final, field_name=var_name)
    except (FileNotFoundError, KeyError) as e:
        print(f"\nERROR: Could not read data: {e}\n")
        sys.exit(1)

    mid_x_index = len(xc) // 2
    mid_y_index = len(yc) // 2
    mid_z_index = len(zc) // 2
    current_axes_row = axes[i]

    # --- Plot 1: Slice at constant X (plotted vs Y) ---
    ax1 = current_axes_row[0]
    slice_init_vs_y = data_init[mid_z_index, :, mid_x_index]
    slice_final_vs_y = data_final[mid_z_index, :, mid_x_index]

    # --- Fit the initial data ---
    guess_init_y = [np.std(slice_init_vs_y), ky_theory, 0, np.mean(slice_init_vs_y)]
    popt_init_y, pcov_init_y = curve_fit(sine_wave_model, yc, slice_init_vs_y, p0=guess_init_y)
    perr_init_y = np.sqrt(np.diag(pcov_init_y))

    # --- Fit the final data ---
    guess_final_y = [np.std(slice_final_vs_y), ky_theory, 0, np.mean(slice_final_vs_y)]
    popt_final_y, pcov_final_y = curve_fit(sine_wave_model, yc, slice_final_vs_y, p0=guess_final_y)
    perr_final_y = np.sqrt(np.diag(pcov_final_y))

    # Plot data and fits
    ax1.plot(yc, slice_init_vs_y, label='Initial Data (t=0)')
    ax1.plot(yc, slice_final_vs_y, label='Final (t= T)', linestyle='dashed')

    # Add text with fit results
    fit_text = (
        f"Initial Fit:\n"
        f"  Amp = {popt_init_y[0]:.4f} ± {perr_init_y[0]:.4f}\n"
        f"  k = {popt_init_y[1]:.4f} ± {perr_init_y[1]:.4f}\n"
        f"Final Fit:\n"
        f"  Amp = {popt_final_y[0]:.4f} ± {perr_final_y[0]:.4f}\n"
        f"  k = {popt_final_y[1]:.4f} ± {perr_final_y[1]:.4f}"
    )
    ax1.text(0.05, 0.35, fit_text, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax1.set_title(f'Slice of {var_name} at x = {xc[mid_x_index]:.3f}')
    ax1.set_xlabel('Y Coordinate')
    ax1.set_ylabel(var_name)
    ax1.legend(loc='lower left')
    ax1.grid(True)

    # --- Plot 2: Slice at constant Y (plotted vs X) ---
    ax2 = current_axes_row[1]
    slice_init_vs_x = data_init[mid_z_index, mid_y_index, :]
    slice_final_vs_x = data_final[mid_z_index, mid_y_index, :]

    # --- Fit the initial data ---
    guess_init_x = [np.std(slice_init_vs_x), kx_theory, 0, np.mean(slice_init_vs_x)]
    popt_init_x, pcov_init_x = curve_fit(sine_wave_model, xc, slice_init_vs_x, p0=guess_init_x)
    perr_init_x = np.sqrt(np.diag(pcov_init_x))

    # --- Fit the final data ---
    guess_final_x = [np.std(slice_final_vs_x), kx_theory, 0, np.mean(slice_final_vs_x)]
    popt_final_x, pcov_final_x = curve_fit(sine_wave_model, xc, slice_final_vs_x, p0=guess_final_x)
    perr_final_x = np.sqrt(np.diag(pcov_final_x))

    # Plot data and fits
    ax2.plot(xc, slice_init_vs_x, label='Initial Data (t=0)')
    ax2.plot(xc, slice_final_vs_x, label='Final (t= T)', linestyle='dashed')

    # Add text with fit results
    fit_text = (
        f"Initial Fit:\n"
        f"  Amp = {popt_init_x[0]:.4f} ± {perr_init_x[0]:.4f}\n"
        f"  k = {popt_init_x[1]:.4f} ± {perr_init_x[1]:.4f}\n"
        f"Final Fit:\n"
        f"  Amp = {popt_final_x[0]:.4f} ± {perr_final_x[0]:.4f}\n"
        f"  k = {popt_final_x[1]:.4f} ± {perr_final_x[1]:.4f}"
    )
    ax2.text(0.75, 0.95, fit_text, transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax2.set_title(f'Slice of {var_name} at y = {yc[mid_y_index]:.3f}')
    ax2.set_xlabel('X Coordinate')
    ax2.legend(loc='lower left')
    ax2.grid(True)

plt.show()

print("\nDone.")
