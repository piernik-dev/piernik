#!/usr/bin/env python3
"""
batch_convert_h5_to_vts.py

A command-line tool that locates all `.h5` files in a given input directory
(and its subdirectories) and invokes a conversion script to produce a `.vts`
file for each one. The output files are saved in a corresponding structure
within the specified output directory.

This wrapper makes no changes to the conversion script; it simply calls it
via a subprocess for each matching input file.

Command-Line Usage
------------------
python batch_convert_h5_to_vts.py <input_directory> <output_directory> <path_to_converter_script>

Example:
python batch_convert_h5_to_vts.py ./raw_data/ ./vtk_data/ ./h5tovts.py
"""

import os
import glob
import subprocess
import sys
import argparse

def find_h5_files(root_dir):
    """Recursively find all .h5 files under root_dir."""
    print(f"Searching for .h5 files in: {root_dir}")
    pattern = os.path.join(root_dir, "**", "*.h5")
    return sorted(glob.glob(pattern, recursive=True))

def make_output_path(in_path, in_root, out_root):
    """
    Given an input path like /in_root/sub/foo.h5,
    return /out_root/sub/foo.vts
    """
    # Create a relative path from the input root
    rel_path = os.path.relpath(in_path, in_root)
    # Split off the extension
    base, _ = os.path.splitext(rel_path)
    # Join with the output root and add the new extension
    return os.path.join(out_root, base + ".vth")

def ensure_dir_exists(path):
    """Create parent directories for a file path if they don’t exist."""
    # Get the directory part of the path
    directory = os.path.dirname(path)
    # If the directory part is not empty and doesn't exist, create it
    if directory and not os.path.isdir(directory):
        os.makedirs(directory, exist_ok=True)
        print(f"Created directory: {directory}")

def main(args):
    """Main execution function driven by command-line arguments."""
    # 1. Locate all HDF5 files
    h5_files = find_h5_files(args.input_dir)
    if not h5_files:
        print(f"Error: No .h5 files found under '{args.input_dir}'", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(h5_files)} files to process.\n")

    # 2. Loop and convert each file
    for idx, h5_path in enumerate(h5_files, start=1):
        # Determine the corresponding output path for the .vts file
        vts_path = make_output_path(h5_path, args.input_dir, args.output_dir)
        
        # Ensure the subdirectory for the output file exists
        ensure_dir_exists(vts_path)

        print(f"[{idx}/{len(h5_files)}] Converting:")
        print(f"    Input:  {h5_path}")
        print(f"    Output: {vts_path}")

        # 3. Invoke the converter script as a subprocess
        try:
            # The command to run: python <converter_script> <input.h5> <output.vts>
            command = [
                sys.executable,        # The current python interpreter
                args.converter_script,
                h5_path,
                vts_path
            ]
            
            result = subprocess.run(
                command,
                check=True,            # Raise an exception if the script fails
                capture_output=True,   # Capture stdout and stderr
                text=True              # Decode stdout/stderr as text
            )
            
            # Print the output from the conversion script for diagnostics
            if result.stdout:
                # Indent the output for clarity
                for line in result.stdout.strip().splitlines():
                    print(f"      > {line}")

        except FileNotFoundError:
            print(f"FATAL ERROR: The converter script was not found at '{args.converter_script}'", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            # This block runs if the converter script returns a non-zero exit code (an error)
            print(f"ERROR: The conversion script failed for '{os.path.basename(h5_path)}'.", file=sys.stderr)
            print("--- Converter's Standard Output ---", file=sys.stderr)
            print(e.stdout, file=sys.stderr)
            print("--- Converter's Standard Error ----", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
            print("-----------------------------------", file=sys.stderr)
        
        print("-" * 20)

    print("Batch conversion complete.")

if __name__ == "__main__":
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(
        description="Batch convert Piernik HDF5 (.h5) files to VTK Structured Grid (.vts) files.",
        formatter_class=argparse.RawTextHelpFormatter # Preserves newlines in help text
    )
    
    parser.add_argument(
        "input_dir",
        help="The root directory containing .h5 files to convert (searched recursively)."
    )
    parser.add_argument(
        "output_dir",
        help="The root directory where the resulting .vts files will be saved."
    )
    parser.add_argument(
        "converter_script",
        help="The path to the Python script that performs the h5-to-vts conversion."
    )

    # Parse the arguments provided by the user
    parsed_args = parser.parse_args()

    # Basic validation of inputs
    if not os.path.isdir(parsed_args.input_dir):
        print(f"Error: Input directory not found at '{parsed_args.input_dir}'", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.isfile(parsed_args.converter_script):
        print(f"Error: Converter script not found at '{parsed_args.converter_script}'", file=sys.stderr)
        sys.exit(1)

    # Run the main function with the parsed arguments
    main(parsed_args)
