#!/usr/bin/env python3
#
# collect positions

# script to parse output of interogate script
# *annotated.tab file and summerise the positions found


import argparse
import os
import csv
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Parse exon data from file", add_help=False)
    file_directory = os.path.realpath(__file__).split("parse_exons.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--file", dest='file',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          default=os.path.join(file_directory, "data", "test.site_proba.csv"),
                          type=str,
                          help="List of input files to be parsed e.g. --file file1.tsv file2.tsv file3.tsv")
 
    optional.add_argument("--output", dest='output',
                          action="store", default="output.tsv",
                          type=str,
                          help="Path to the output file (default: output.tsv)")
    
    return parser.parse_args()

def parse_file(file_path, data):
    """
    Parse exon data from a given file and update the provided dictionary.

    This function reads a tab-separated values (TSV) file containing exon data,
    and updates the given defaultdict with transcript IDs, exon numbers, and positions.

    Args:
        file_path (str): The path to the input file to be parsed.
        data (defaultdict): A defaultdict where each key is a transcript ID,
                            and each value is another defaultdict. The inner
                            defaultdict's keys are exon numbers, and the values
                            are lists of positions.

    Returns:
        None: The function updates the provided defaultdict in place.
    """
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        # Dynamically determine if 'position' or 'positions' is in the file
        # for weird reasons I was getting both names. Maybe an old scripts output?
        column_name = 'position' if 'position' in csv_reader.fieldnames else 'positions'

        for row in csv_reader:
            transcript_id = row['transcript_id']
            exon_number = row['exon_number']
            position = row[column_name]  # Use the dynamically determined name
            data[transcript_id][exon_number].append(position)

def main():
    args = get_args()
    file_paths = args.file
    output_file = args.output

    # WARNING: adds to the same dictionary for all the input files. 
    
    all_data = defaultdict(lambda: defaultdict(list))
    
    for file_path in file_paths:
        parse_file(file_path, all_data)
    
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['transcript_id', 'exon_number', 'positions', 'num_positions', 
                         'num_unique_positions', 'summary'])
        
        for transcript_id, exons in all_data.items():
            for exon_number, positions in exons.items():
                num_positions = len(positions)
                unique_positions = sorted(set(positions), key=int)
                num_unique_positions = len(unique_positions)
                positions.sort(key=int)
                positions_str = ','.join(positions)
                summary = ', '.join([f"{pos}({positions.count(pos)})" for pos in unique_positions])
                writer.writerow([transcript_id, exon_number, positions_str, num_positions, num_unique_positions, summary])
    
    print(f"Data successfully written to {output_file}")

if __name__ == "__main__":
    main()

