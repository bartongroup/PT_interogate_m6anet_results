import argparse
import os
import csv
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(description="Compare exon modification data from multiple files", add_help=False)
    file_directory = os.path.realpath(__file__).split("compare_exons.py")[0]
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument("--files", dest='files',
                          action="store",
                          nargs='+',  # This allows multiple arguments
                          required=True,
                          default=[os.path.join(file_directory, "data", "test1.tsv"), os.path.join(file_directory, "data", "test2.tsv")],
                          type=str,
                          help="List of input files to be parsed and compared e.g. --files file1.tsv file2.tsv")
 
    optional.add_argument("--output", dest='output',
                          action="store", default=None,
                          type=str,
                          help="Path to the output file (default: derived from input filenames)")
    
    return parser.parse_args()

def parse_file(file_path):
    """
    Parse exon data from a given file and return a dictionary.

    This function reads a tab-separated values (TSV) file containing exon data,
    and returns a dictionary with transcript IDs, exon numbers, and positions.

    Args:
        file_path (str): The path to the input file to be parsed.

    Returns:
        dict: A dictionary where each key is a transcript ID, and each value
              is another dictionary. The inner dictionary's keys are exon
              numbers, and the values are sets of unique positions.
    """
    data = defaultdict(lambda: defaultdict(set))  # Using set to ensure unique positions
    
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            transcript_id = row['transcript_id']
            exon_number = row['exon_number']
            positions = row['positions'].split(',')
            positions = [pos.strip() for pos in positions if pos.strip()]  # Clean and filter positions
            data[transcript_id][exon_number].update(positions)
            print(f"Parsed {transcript_id} {exon_number} with positions: {positions}")  # Debug print
    
    print(f"Data from {file_path}: {dict(data)}")  # Debug print entire data
    return data

def compare_files(data1, data2, filename1, filename2):
    """
    Compare exon modification data from two files.

    This function takes two dictionaries containing exon modification data
    and compares the transcripts and positions of modifications.

    Args:
        data1 (dict): A dictionary containing exon modification data from the first file.
        data2 (dict): A dictionary containing exon modification data from the second file.
        filename1 (str): The name of the first file.
        filename2 (str): The name of the second file.

    Returns:
        dict: A dictionary containing the comparison results.
    """
    comparison_results = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    
    all_transcripts = set(data1.keys()).union(data2.keys())
    
    for transcript_id in all_transcripts:
        all_exons = set(data1[transcript_id].keys()).union(data2[transcript_id].keys())
        
        for exon_number in all_exons:
            positions1 = data1[transcript_id][exon_number] if exon_number in data1[transcript_id] else set()
            positions2 = data2[transcript_id][exon_number] if exon_number in data2[transcript_id] else set()
            
            common_positions = positions1.intersection(positions2)
            unique_positions1 = positions1 - common_positions
            unique_positions2 = positions2 - common_positions

            comparison_results[transcript_id][exon_number]['common'] = sorted(list(common_positions), key=int)
            comparison_results[transcript_id][exon_number][f'unique_{os.path.basename(filename1)}'] = sorted(list(unique_positions1), key=int)
            comparison_results[transcript_id][exon_number][f'unique_{os.path.basename(filename2)}'] = sorted(list(unique_positions2), key=int)

            print(f"Transcript {transcript_id} Exon {exon_number}: Common positions: {common_positions}")  # Debug print
            print(f"Unique positions for {transcript_id} {exon_number} in file {filename1}: {unique_positions1}")  # Debug print
            print(f"Unique positions for {transcript_id} {exon_number} in file {filename2}: {unique_positions2}")  # Debug print
    
    return comparison_results

def main():
    args = get_args()
    file_paths = args.files
    output_file = args.output
    
    if output_file is None:
        base_names = [os.path.splitext(os.path.basename(f))[0] for f in file_paths]
        output_file = "_vs_".join(base_names) + "_comparison.tsv"
    
    # Parse files separately
    data1 = parse_file(file_paths[0])
    data2 = parse_file(file_paths[1])
    
    # Compare data from the two files
    comparison_results = compare_files(data1, data2, file_paths[0], file_paths[1])
    
    # Write comparison results to the output file
    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        header = ['transcript_id', 'exon_number', 'common_positions', 'num_common_positions',
                  f'unique_positions_{os.path.basename(file_paths[0])}', f'num_unique_positions_{os.path.basename(file_paths[0])}',
                  f'unique_positions_{os.path.basename(file_paths[1])}', f'num_unique_positions_{os.path.basename(file_paths[1])}']
        writer.writerow(header)
        
        for transcript_id, exons in comparison_results.items():
            for exon_number, positions_data in exons.items():
                common_positions = positions_data['common']
                num_common_positions = len(common_positions)
                row = [transcript_id, exon_number, ', '.join(common_positions), num_common_positions]
        
                unique_positions1 = positions_data.get(f'unique_{os.path.basename(file_paths[0])}', [])
                num_unique_positions1 = len(unique_positions1) if unique_positions1 else 0
                unique_positions2 = positions_data.get(f'unique_{os.path.basename(file_paths[1])}', [])
                num_unique_positions2 = len(unique_positions2) if unique_positions2 else 0
        
                unique_positions1_str = ', '.join(unique_positions1) if unique_positions1 else '0'
                unique_positions2_str = ', '.join(unique_positions2) if unique_positions2 else '0'
        
                row.extend([unique_positions1_str, num_unique_positions1,
                            unique_positions2_str, num_unique_positions2])
        
                writer.writerow(row)
    
    print(f"Comparison data successfully written to {output_file}")

if __name__ == "__main__":
    main()
