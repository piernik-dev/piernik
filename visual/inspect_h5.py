import h5py
import sys

def print_hdf5_structure(name, obj):
    """Callback function to print the structure of an HDF5 file."""
    indent = '  ' * name.count('/')
    print(f"{indent}--> {name.split('/')[-1]}")
    # Also print attributes if any
    if obj.attrs:
        for key, val in obj.attrs.items():
            print(f"{indent}    - ATTR: {key} = {val}")

def inspect_hdf5_file(filename):
    """
    Prints the full group, dataset, and attribute structure of an HDF5 file.
    """
    print("-" * 60)
    print(f"Inspecting structure of file: {filename}")
    print("-" * 60)
    try:
        with h5py.File(filename, 'r') as f:
            print(f"ROOT GROUP ('/')")
            # Print attributes of the root group first
            if f.attrs:
                for key, val in f.attrs.items():
                    print(f"    - ATTR: {key} = {val}")
            
            # Walk through all other objects in the file
            f.visititems(print_hdf5_structure)
    except Exception as e:
        print(f"Could not read file. Error: {e}")
    print("-" * 60)

# --- Main part of the script ---
if __name__ == "__main__":
    # If a filename is given on the command line, use it.
    # Otherwise, default to '0000.h5'.
    if len(sys.argv) > 1:
        file_to_inspect = sys.argv[1]
    else:
        file_to_inspect = '0000.h5' 
        
    inspect_hdf5_file(file_to_inspect)
