#!/usr/bin/env python3

import sys
import subprocess
import os

def get_script_path(script_name):
    package_root = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(package_root, 'iggypop', script_name)

def main():
    # Check if at least one argument is provided (the run type)
    if len(sys.argv) < 2:
        print("Error: Run type ('cds', 'gb', 'format', 'gagga', or 'primers') not specified.")
        sys.exit(1)

    # The first positional argument is the run type
    runtype = sys.argv[1]

    # Verify the run type is correct
    if runtype not in ['cds', 'gb', 'format', 'gagga', 'primers']:
        print("Error: Invalid run type specified. Choose 'cds', 'gb', 'format', 'gagga', or 'primers'.")
        sys.exit(1)

    # Map runtype to script names
    script_map = {
        'cds': 'iggypop_cds.py',
        'gb': 'iggypop_gb.py',
        'format': 'iggypop_format.py',
        'gagga': 'iggypop_gagga.py',
        'primers': 'indexing_primers.py'
    }

    script_to_run = script_map[runtype]
    script_path = get_script_path(script_to_run)

    # Additional arguments to be passed can be collected from sys.argv[2:]
    additional_args = sys.argv[2:]

    # If no additional arguments are provided, use --help
    if not additional_args:
        additional_args = ['--help']

    # Building the command to run the subprocess
    command = ['python', script_path] + additional_args

    # Execute the subprocess with the selected script and additional arguments
    subprocess.run(command)

if __name__ == "__main__":
    main()
