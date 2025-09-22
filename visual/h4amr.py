#!/usr/bin/env python3
import os
import sys
import struct
import h5py
import numpy as np

# ----------------------------
# Helpers to read block attrs
# ----------------------------

def _get_attr(attrs, *names, required=False):
    """Try a list of attribute names and return the first that exists."""
    for n in names:
        if n in attrs:
            return attrs[n]
    if required:
        raise KeyError(f"Required attribute not found. Tried: {names}")
    return None

def _normalize_fbnd(fbnd):
    """fbnd may be stored as shape (3,2) or (2,3); return (lo[3], hi[3])."""
    a = np.asarray(fbnd)
    if a.shape == (3, 2):
        lo = a[:, 0]
        hi = a[:, 1]
    elif a.shape == (2, 3):
        lo = a[0, :]
        hi = a[1, :]
    else:
        raise ValueError(f"fbnd has unexpected shape {a.shape}, expected (3,2) or (2,3)")
    return np.array(lo, dtype=np.float64), np.array(hi, dtype=np.float64)

def _block_geo(attrs, global_origin=None):
    """
    Return per-block (dims, origin, spacing) from block attributes.
    dims are cell counts (nx,ny,nz).
    origin and spacing are physical (double).
    """
    # dims (cell counts)
    dims = _get_attr(attrs, "n_b", "nb", "dims", required=True)
    dims = np.array(dims, dtype=int)
    if dims.size != 3:
        raise ValueError(f"n_b/nb/dims has wrong length: {dims}")

    # Prefer fbnd if present
    fbnd = _get_attr(attrs, "fbnd", "bnd", "bounds", required=False)
    if fbnd is not None:
        lo, hi = _normalize_fbnd(fbnd)
        spacing = (hi - lo) / np.maximum(1, dims.astype(np.float64))
        origin = lo
        return dims, origin, spacing

    # Else use dl (spacing) + off (cell offset) relative to a global origin
    dl = _get_attr(attrs, "dl", "dcell", "spacing", required=False)
    off = _get_attr(attrs, "off", "offset", required=False)
    if dl is None or off is None:
        raise KeyError("Neither fbnd nor (dl+off) found in block attributes")

    dl = np.array(dl, dtype=np.float64)
    off = np.array(off, dtype=int)
    if global_origin is None:
        raise KeyError("Global origin is required to compute block origin from off*dl")
    origin = np.array(global_origin, dtype=np.float64) + off.astype(np.float64) * dl
    spacing = dl
    return dims, origin, spacing

def _global_origin_from_root(f):
    # Most PIERNIK outputs provide domain_left_edge under simulation parameters
    grp = _get_attr(f, "simulation_parameters")
    if grp is None:
        # Try ‘/’ attrs
        root_attrs = f.attrs
        return _get_attr(root_attrs, "domain_left_edge", "xmin", "domain_origin", required=True)
    else:
        return _get_attr(f["simulation_parameters"].attrs, "domain_left_edge", "xmin", "domain_origin", required=True)

# ----------------------------
# VTI / VTM writers
# ----------------------------

def write_vti(path, dims, origin, spacing, arrays):
    """
    Write a single .vti (XML ImageData with appended raw) for one AMR block.
    arrays: dict name -> np.ndarray with shape (nz,ny,nx), dtype float32/float64.
    Stored as CellData (VTK expects UInt32 length headers, little-endian float32).
    """
    nx, ny, nz = map(int, dims)  # dims are cell counts
    # VTK Extent is in *points*: to hold nx*ny*nz CELLS, we must use 0..nx, 0..ny, 0..nz
    whole_extent = (0, nx, 0, ny, 0, nz)
    piece_extent = whole_extent

    # Only keep arrays with correct shape
    kept = {}
    for name, arr in arrays.items():
        A = np.asarray(arr)
        if A.ndim != 3:
            # skip vector/tensor here; (optional) expand into components upstream if needed
            continue
        if A.shape != (nz, ny, nx):
            raise ValueError(f"Array '{name}' has shape {A.shape}, expected (nz,ny,nx)=({nz},{ny},{nx}).")
        kept[name] = A.astype(np.float32, copy=False)

    names = list(kept.keys())
    n_arrays = len(names)
    if n_arrays == 0:
        return  # nothing to write

    # Each raw block is len = nx*ny*nz*4 bytes (Float32), preceded by 4-byte UInt32
    payload_len = nx * ny * nz * 4
    offsets = [0]
    for i in range(1, n_arrays):
        offsets.append(offsets[-1] + 4 + payload_len)

    # 1) XML header (formatted)
    with open(path, "wb") as fxml:
        fxml.write(b'<?xml version="1.0"?>\n')
        fxml.write(b'<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt32">\n')
        fxml.write(
            f'  <ImageData WholeExtent="{whole_extent[0]} {whole_extent[1]} {whole_extent[2]} {whole_extent[3]} {whole_extent[4]} {whole_extent[5]}" '
            f'Origin="{origin[0]:.16e} {origin[1]:.16e} {origin[2]:.16e}" '
            f'Spacing="{spacing[0]:.16e} {spacing[1]:.16e} {spacing[2]:.16e}">\n'.encode()
        )
        fxml.write(
            f'    <Piece Extent="{piece_extent[0]} {piece_extent[1]} {piece_extent[2]} {piece_extent[3]} {piece_extent[4]} {piece_extent[5]}">\n'.encode()
        )
        fxml.write(b'      <PointData/>\n')
        fxml.write(b'      <CellData>\n')
        for name, off in zip(names, offsets):
            # Use names as-is; if you want GDF-style renaming, map here
            line = f'        <DataArray type="Float32" Name="{name}" format="appended" offset="{off}"/>\n'
            fxml.write(line.encode())
        fxml.write(b'      </CellData>\n')
        fxml.write(b'    </Piece>\n')
        fxml.write(b'  </ImageData>\n')
        fxml.write(b'  <AppendedData encoding="raw">\n')
        fxml.write(b'_')  # appended-data marker

        # 2) Appended RAW (UInt32 length + little-endian Float32 data, x-fastest)
        for name in names:
            fxml.write(struct.pack("<I", payload_len))  # 4-byte length
            arr = kept[name]
            # Our arrays are (nz,ny,nx) with x fastest in C order -> OK for VTK (x,y,z)
            fxml.write(arr.tobytes(order="C"))

        # 3) Footer
        fxml.write(b'\n  </AppendedData>\n')
        fxml.write(b'</VTKFile>\n')

def write_vtm(path, vti_filenames):
    with open(path, "wb") as f:
        f.write(b'<?xml version="1.0"?>\n')
        f.write(b'<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">\n')
        f.write(b'  <vtkMultiBlockDataSet>\n')
        for ib, vti in enumerate(vti_filenames):
            line = f'    <DataSet index="{ib}" name="block_{ib}" file="{os.path.basename(vti)}"/>\n'
            f.write(line.encode())
        f.write(b'  </vtkMultiBlockDataSet>\n')
        f.write(b'</VTKFile>\n')

# ----------------------------
# H5 -> VTM (with .vti blocks)
# ----------------------------

def convert_h5_to_vtm(h5_path, fields=None, out_dir=None):
    """
    Convert a PIERNIK HDF5 output (possibly AMR) into a .vtm with one .vti per block.
    If 'fields' is None, include every dataset inside each block group.
    """
    if out_dir is None:
        out_dir = os.path.dirname(os.path.abspath(h5_path))
    base = os.path.splitext(os.path.basename(h5_path))[0]
    vtm_path = os.path.join(out_dir, base + ".vtm")

    vti_paths = []
    with h5py.File(h5_path, "r") as f:
        data_grp = f["data"]
        block_names = sorted(list(data_grp.keys()))

        # Global origin (only needed for the off*dl path)
        try:
            global_origin = _global_origin_from_root(f)
        except Exception:
            global_origin = None  # okay if every block has fbnd

        for ib, bname in enumerate(block_names):
            g = data_grp[bname]
            attrs = g.attrs

            # Decide which arrays to include for this block
            if fields is None:
                # Take all datasets directly in this block group that look scalar and 3D
                block_fields = [k for k, ds in g.items() if isinstance(ds, h5py.Dataset) and ds.ndim == 3]
            else:
                block_fields = [k for k in fields if k in g]

            if not block_fields:
                # skip empty block
                continue

            dims, origin, spacing = _block_geo(attrs, global_origin=global_origin)
            nx, ny, nz = map(int, dims)
            arrays = {}
            for k in block_fields:
                A = np.array(g[k], dtype=np.float32)
                # Expect (z,y,x) order from your writer
                if A.shape != (nz, ny, nx):
                    raise ValueError(f"Block {bname}: dataset '{k}' shape {A.shape} != expected (nz,ny,nx)=({nz},{ny},{nx})")
                arrays[k] = A

            vti_path = os.path.join(out_dir, f"{base}_block_{ib}.vti")
            write_vti(vti_path, dims, origin, spacing, arrays)
            vti_paths.append(vti_path)

    # Write the multiblock
    write_vtm(vtm_path, vti_paths)
    return vtm_path, vti_paths

# ----------------------------
# CLI
# ----------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python h5_to_vtm_amr.py <file.h5> [field1 field2 ...]")
        sys.exit(1)
    h5_file = sys.argv[1]
    fields = sys.argv[2:] if len(sys.argv) > 2 else None
    vtm, vtis = convert_h5_to_vtm(h5_file, fields=fields)
    print("Wrote:", vtm)
    for p in vtis:
        print("  ", p)

