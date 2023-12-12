#!/usr/bin/env python3

# Author(s): Roderick S.K. Westerman & Sander W. van der Laan | github.com/swvanderlaan | s.w.vanderlaan@gmail.com
# Date: 2023-12-11
# License: CC-BY-NC-ND 4.0 International
# Disclaimer: "The script is given AS-IS. No liability is taken for its use."
# Version: 1.0
# Usage: ./loci_lookup.py --fileA <path_to_file_A> --fileB <path_to_file_B> --columnA Chr BP --columnB chromosome base_pair_location [--log <path_to_log_file>] [--verbose] 
# Description: This script takes two input files (file A and file B) and performs a lookup based on specified columns.
#              The result is written to an output file, and information about the process is logged.
#              The script dynamically installs required packages if they are not installed, checks if the input files are gzipped,
#              and detects the delimiter used in the input files.

"""
Loci Lookup Script

This script takes two input files (file A and file B) and performs a lookup based on specified columns.
The result is written to an output file, and information about the process is logged.

Usage:
  ./loci_lookup.py --fileA <path_to_file_A> --fileB <path_to_file_B> --columnA Chr BP --columnB chromosome base_pair_location [--log <path_to_log_file>] [--verbose] 

Options:
  -A <path_to_file_A>, --fileA <path_to_file_A>           Path to file A (non-gzipped; space or tab-separated).
  -B <path_to_file_B>, --fileB <path_to_file_B>           Path to file B (gzipped; space or tab-separated).
  -cA CHR POS, --columnA CHR POS                           Columns in file A.
  -cB chromosome basepair_position, --columnB chromosome basepair_position  Columns in file B.
  -dA s --delimiterA s                                    Delimiter for file A (s[pace], t[ab], c[omma], [se]m[icolon], [co]l[on]). [NOT IMPLEMENTED YET]
  -db s --delimiterB s                                    Delimiter for file B (s[pace], t[ab], c[omma], [se]m[icolon], [co]l[on]). [NOT IMPLEMENTED YET]

  -o <output_directory>, --output <output_directory>     Output directory (default: current working directory).
  -l <path_to_log_file>, --log <path_to_log_file>        Path to the log file (default: OUTPUT/loci_lookup.log).
  -v, --verbose                                          Enable verbose mode.
  -V, --version                                          Show version information.
  -h, --help                                             Show this help message and exit.
"""

import argparse
import os
import logging
# gzip is also imported if file ends in .gz

# Check if the required packages are installed
try:
    import pandas as pd
except ModuleNotFoundError:
    ### This part is a security concern
    # import subprocess     
    # subprocess.run(["pip", "install", 'pandas'])
    # import pandas as pd
    ### This is an alternative.
    print("Please install pandas 'pip install pandas'")
    raise

VERSION = '1.0'
COPYRIGHT = 'Copyright 1979-2023 | CC-BY-NC-ND License | Created by: Roderick S.K. Westerman & Sander W. van der Laan | github.com/swvanderlaan'
DEFAULT_OUTPUT_FILE = 'loci_lookup.txt'

# Setup logging
def setup_logger(log_file):
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# Detect delimiter
def detect_delimiter(file_path, num_lines=5):
    """num_lines to read a couple of lines to find the delimiter"""
    file_ext = os.path.splitext(file_path)[1]
    potential_delimiters = [',', '\t', ';']  # Add more if needed

    if file_ext != ".gz":
        with open(file_path, 'r') as file:
            lines = [file.readline().strip() for _ in range(num_lines)]
    else:
        import gzip
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            lines = [file.readline().strip() for _ in range(num_lines)]

    for delimiter in potential_delimiters:
        if all(delimiter in line for line in lines):
            return delimiter.encode('unicode_escape').decode('utf-8')
    return None

def main():
    print(f'============================')
    print(f'     Loci Lookup Script')
    print(f'        v{VERSION}')
    print(f'============================')
    print(f'')
    print(f'Beginning lookup process...')
    
    parser = argparse.ArgumentParser(description=f'Loci Lookup Script v{VERSION}. Lookup variants in file B based on file A. The lookup is based on specified columns designating the chromosome and basepair position.')

    parser.add_argument('-A', '--fileA', required=True, help='Path to non-gzipped file A.')
    parser.add_argument('-B', '--fileB', required=True, help='Path to gzipped file B.')
    parser.add_argument('-cA', '--columnA', required=True, nargs=2, metavar=('CHR','POS'), help='Columns in file A (space or tab-separated).')
    parser.add_argument('-cB', '--columnB', required=True, nargs=2, metavar=('chromosome', 'basepair_position'), help='Columns in file B (space or tab-separated).')
    parser.add_argument('-dA', '--delimiterA', required=False, default='t', help='Delimiter for file A (s[pace], t[ab], c[omma], [se]m[icolon], [co]l[on]). [NOT IMPLEMENTED YET]')
    parser.add_argument('-dB', '--delimiterB', required=False, default='t', help='Delimiter for file B (s[pace], t[ab], c[omma], [se]m[icolon], [co]l[on]). [NOT IMPLEMENTED YET]')
    parser.add_argument('-o', '--output', default='', help='Output directory + filename (default: current working directory).')
    parser.add_argument('-l', '--log', required=False, default='./loci_lookup.log', help='Path to the log file (default: loci_lookup.log).')
    parser.add_argument('-p', '--progress', action='store_true', help='Show progress bar.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode.', default=False)
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s v{VERSION}')

    args = parser.parse_args()
    
    # Set the output directory
    current_wd = os.getcwd() # Current working directory
    if args.verbose:
        print(f'- Current working directory: [{current_wd}]')
        logging.info(f'Current working directory: {current_wd}')

    # Set up the logging
    # Get the filename and directory from the log argument
    log_filename = os.path.basename(args.log)
    log_dir_path = os.path.dirname(args.log)
    # Create LOG directory if it doesn't exist
    log_dir = os.path.join(current_wd, log_dir_path)
    os.makedirs(os.path.join(log_dir_path), exist_ok=True)
        
    # Verbose: print the log directory
    if args.verbose:
        print(f'- Created log directory (if it didn\'t exist): [{log_dir}]')
        print(f'- Saving log of lookup: [{log_filename}]')
        logging.info(f'Creating output directory: {log_dir}')
        logging.info(f'Saving log of lookup: {log_filename}')

    # Set up the file path for the log file
    log_filename_path = os.path.join(log_dir, log_filename) if args.log else os.path.join(current_wd, log_filename)

    # Verbose: print the log file path
    if args.verbose:
        print(f'- Log file: [{log_filename_path}]')
        logging.info(f'Log file: {log_filename_path}')

    setup_logger(log_filename_path)

    # Set the output filename and directory
    if args.output:
        # get the filename and directory from the output argument
        output_filename = os.path.basename(args.output)
        output_dir_path = os.path.dirname(args.output)
        # Create OUTPUT directory if it doesn't exist
        output_dir = os.path.join(current_wd, output_dir_path) # os.path.abspath(os.path.join(current_wd, output_dir_path))
        os.makedirs(os.path.join(output_dir), exist_ok=True)
        
        # Verbose: print the output directory
        if args.verbose:
            print(f'- Created output directory (if it didn\'t exist): [{output_dir}]')
            print(f'- Saving looked up data in: [{output_filename}]')
            logging.info(f'Creating output directory: {output_dir}')
            logging.info(f'Saving looked up data in: {output_filename}')

    # Verbose: print the input files
    if args.verbose:
        print(f'- Loci of interest in File A: [{args.fileA}]')
        print(f'- Look-up data in File B: [{args.fileB}]')
        logging.info(f'Loci of interest in File A: {args.fileA}')
        logging.info(f'Look-up data in File B: {args.fileB}')
    
    # Detect delimiters for non-gzipped files (can be adjusted if necessary)
    # delimiter_A = '\t'  # Assuming fileA is tab-delimited
    # delimiter_B = '\t'  # Assuming fileB is tab-delimited    
    delimiter_A = detect_delimiter(args.fileA) # Slow
    delimiter_B = detect_delimiter(args.fileB) # Slow
    
    # Verbose: print the detected delimiters
    if args.verbose:
        detect_delimiter_A = detect_delimiter(args.fileA)
        print(f"- Delimiter for file A: [{detect_delimiter_A}].")
        detect_delimiter_B = detect_delimiter(args.fileB)
        print(f"- Delimiter for file B: [{detect_delimiter_B}].")
        logging.info(f"Delimiter for file B: [{detect_delimiter_A}].")
        logging.info(f"Delimiter for file B: [{detect_delimiter_B}].")

    # Set the output file path
    output_file_path = os.path.join(output_dir, output_filename) if args.output else os.path.join(current_wd, DEFAULT_OUTPUT_FILE)

    # Verbose: print the output file path
    if args.verbose:
        print(f'- Output file: [{output_file_path}]')
        logging.info(f'Output file: {output_file_path}')

    # Read the files    
    df_fileA = pd.read_csv(args.fileA, sep=delimiter_A, engine='python')
    df_fileB = pd.read_csv(args.fileB, sep=delimiter_B, engine='python')
    
    # Verbose: print the headers of the files
    if args.verbose:
        print(f'- File A header: {df_fileA.columns}')
        print(f'- File B header: {df_fileB.columns}')
        logging.info(f'File A header: {df_fileA.columns}')
        logging.info(f'File B header: {df_fileB.columns}')
    
    # Merge the files
    if args.verbose:
        # Alternatively, write the merged file to disk with all columns
        print(f'- Saving all columns of the merged file to disk while keeping all rows from fileA (useful when debugging in case of \'missing\' variants).')
        # Actually merge the files
        df_merged = pd.merge(df_fileA, df_fileB, left_on=args.columnA, right_on=args.columnB, how='left')
        # Write the merged file to disk
        df_merged.to_csv(output_file_path, sep='\t', index=False)
    else:
        # Actually merge the files
        df_merged = pd.merge(df_fileA, df_fileB, left_on=args.columnA, right_on=args.columnB, how='inner')
        # Write the merged file to disk
        df_merged.to_csv(output_file_path, sep='\t', index=False, columns=df_fileB.columns)

    # Verbose: print the contents of the output file
    if args.verbose:
        print(f'- Checking contents of output file...')
        print(df_merged.head())
        logging.info(f'Checking contents of output file...')
        logging.info(df_merged.head())
        pass
    
    return df_merged
    
if __name__ == "__main__":
    from time import time
    from datetime import timedelta
    # Start the timer
    t1 = time()
    # Run the main function
    df = main()
    # Stop the timer
    t2 = time()
    # Print the time taken
    print(f'Time taken: {timedelta(seconds=t2-t1)}')
    print(f'')
    print(f'============================')
    print(f'{COPYRIGHT}')