
import os
import pandas as pd
import argparse
import logging
from scipy.stats import fisher_exact


# Define thresholds as constants
INDIV_THRESHOLD = 0.0032978046219796  # For data.indiv_proba.csv
SITE_THRESHOLD = 0.9  # For data.site_proba.csv

# Set up logging with both terminal and file output
log_filename = "intron_retention_m6anet_methylated.log"

# Create a logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)  # Set to INFO or DEBUG as needed

# Create handlers
console_handler = logging.StreamHandler()  # For terminal output
file_handler = logging.FileHandler(log_filename)  # For file output

# Set levels for handlers if different levels are needed (optional)
console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

# Add handlers to the logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

def load_and_filter_main_data(main_data_path):
    """
    Load and filter the main data file for rows with 'intron_retention' in the 'subcategory' column.
    
    Parameters:
    main_data_path (str): Path to the main data file containing isoform, gene, and subcategory information.

    Returns:
    pd.DataFrame: A DataFrame containing filtered data with columns for isoform, associated_gene, and associated_transcript.
    """
    main_df = pd.read_csv(main_data_path, sep='\t')
    filtered_df = main_df[main_df['subcategory'] == 'intron_retention'][['isoform', 
                                                                         'associated_gene', 
                                                                         'associated_transcript']]
    logger.info(f"Filtered main data for intron_retention. Retained {filtered_df.shape[0]} rows.")
    return filtered_df


def filter_m6anet_files(m6anet_dir, mapped_data, match_criteria="original"):
    """
    Filter m6anet files based on specified thresholds for probability values,
    allowing selection between original, strict, and gene-strict matching criteria.

    Parameters:
    m6anet_dir (str): Path to the directory containing subdirectories with m6anet output files.
    mapped_data (pd.DataFrame): A DataFrame with isoform, associated_gene, associated_transcript, and read_index.
    match_criteria (str): The criteria to use for filtering. Options are "original", "strict", or "gene_strict".

    Returns:
    tuple: DataFrames for significant and non-significant entries based on the selected matching criteria.
    """
    all_indiv_significant = []
    all_indiv_non_significant = []
    total_non_sig_rows = 0

    for subdir in os.listdir(m6anet_dir):
        subdir_path = os.path.join(m6anet_dir, subdir)
        if os.path.isdir(subdir_path):
            indiv_path = os.path.join(subdir_path, 'data.indiv_proba.csv')

            if os.path.exists(indiv_path):
                indiv_df = pd.read_csv(indiv_path)
                
                # Split data into significant and non-significant based on threshold
                indiv_significant = indiv_df[indiv_df['probability_modified'] < INDIV_THRESHOLD]
                indiv_non_significant = indiv_df[indiv_df['probability_modified'] >= INDIV_THRESHOLD]

                # Log the current merge attempt
                logger.info(f"Attempting merge for {subdir} with criteria {match_criteria}")
                logger.info(f"Columns in mapped_data: {mapped_data.columns}")
                logger.info(f"Columns in indiv_significant: {indiv_significant.columns}")

                # Perform merge based on matching criteria
                if match_criteria == "strict":
                    # Strict transcript matching (read_index + transcript_id == associated_transcript)
                    if 'transcript_id' in indiv_significant.columns:
                        orig_sig = mapped_data.merge(indiv_significant, left_on=['read_index', 'associated_transcript'],
                                                     right_on=['read_index', 'transcript_id'], how='inner')
                        orig_non_sig = mapped_data.merge(indiv_non_significant, left_on=['read_index', 'associated_transcript'],
                                                         right_on=['read_index', 'transcript_id'], how='inner')
                    else:
                        logger.warning(f"transcript_id column missing in {subdir}, skipping strict merge.")
                        continue
                elif match_criteria == "gene_strict":
                    # Gene strict matching (read_index + associated_gene)
                    if 'associated_gene' in mapped_data.columns:
                        orig_sig = mapped_data.merge(indiv_significant, on=['read_index'], how='inner')
                        orig_sig = orig_sig[orig_sig['associated_gene'] == orig_sig['transcript_id'].str.split('.').str[0]]
                        
                        orig_non_sig = mapped_data.merge(indiv_non_significant, on=['read_index'], how='inner')
                        orig_non_sig = orig_non_sig[orig_non_sig['associated_gene'] == orig_non_sig['transcript_id'].str.split('.').str[0]]
                    else:
                        logger.warning(f"associated_gene column missing in mapped data, skipping gene strict merge.")
                        continue
                else:  # Original matching (relaxed matching allowing novel IDs)
                    orig_sig = mapped_data.merge(indiv_significant, on='read_index', how='inner')
                    orig_sig = orig_sig[
                        (orig_sig['associated_transcript'] == orig_sig['transcript_id']) |
                        ((orig_sig['associated_transcript'] == "novel") & 
                         (orig_sig['associated_gene'] == orig_sig['transcript_id'].str.split('.').str[0]))
                    ]
                    orig_non_sig = mapped_data.merge(indiv_non_significant, on='read_index', how='inner')
                    orig_non_sig = orig_non_sig[
                        (orig_non_sig['associated_transcript'] == orig_non_sig['transcript_id']) |
                        ((orig_non_sig['associated_transcript'] == "novel") & 
                         (orig_non_sig['associated_gene'] == orig_non_sig['transcript_id'].str.split('.').str[0]))
                    ]

                # Log and accumulate results
                if not orig_non_sig.empty:
                    logger.info(f"Non-significant data retained for {subdir} - {orig_non_sig.shape[0]} rows")
                    total_non_sig_rows += orig_non_sig.shape[0]  # Add to the total count
                else:
                    logger.warning(f"No non-significant data found for {subdir}")

                all_indiv_significant.append(orig_sig)
                all_indiv_non_significant.append(orig_non_sig)

    # Log the cumulative total of non-significant rows across all files
    logger.info(f"Total non-significant rows across all files: {total_non_sig_rows}")

    # Combine all filtered data for each output
    combined_significant = pd.concat(all_indiv_significant, ignore_index=True) if all_indiv_significant else pd.DataFrame()
    combined_non_significant = pd.concat(all_indiv_non_significant, ignore_index=True) if all_indiv_non_significant else pd.DataFrame()

    return combined_significant, combined_non_significant




def map_read_id_to_index(filtered_data, index_to_reads_dir):
    """
    Map isoform (read ID) to read_index using files in the index_to_reads directory.
    
    Parameters:
    filtered_data (pd.DataFrame): A DataFrame containing filtered isoform data.
    index_to_reads_dir (str): Path to the directory containing index to reads mapping files.

    Returns:
    pd.DataFrame: A DataFrame containing isoform, associated_gene, associated_transcript, and read_index.
    """
    all_index_dfs = []
    
    for filename in os.listdir(index_to_reads_dir):
        if filename.endswith(".tsv"):
            file_path = os.path.join(index_to_reads_dir, filename)
            try:
                index_df = pd.read_csv(file_path, sep='\t')
                if index_df.empty:
                    logger.warning(f"{filename} is empty and will be skipped.")
                    continue
                index_df = index_df[['read_name', 'read_index']].rename(columns={'read_name': 'isoform'})
                merged_df = filtered_data.merge(index_df, on='isoform', how='inner')
                all_index_dfs.append(merged_df)
                logger.info(f"Processed {filename} for read index mapping. Merged rows: {merged_df.shape[0]}.")
            except pd.errors.EmptyDataError:
                logger.warning(f"{filename} is empty or has no columns to parse. Skipping this file.")
    
    combined_df = pd.concat(all_index_dfs, ignore_index=True)
    logger.info(f"Combined all mapped data frames. Total rows after merging: {combined_df.shape[0]}.")
    return combined_df


def calculate_summary_stats(filtered_main, indiv_significant, indiv_non_significant, output_file):
    """
    Calculate and write summary statistics for significant and non-significant data to an output file.
    
    Parameters:
    filtered_main (pd.DataFrame): DataFrame containing filtered data from main file.
    indiv_significant (pd.DataFrame): DataFrame with data passing the threshold.
    indiv_non_significant (pd.DataFrame): DataFrame with data failing the threshold.
    output_file (file handle): File handle for saving the summary output.

    Returns:
    None
    """
    # Total transcripts with retained introns
    total_retained_introns = len(filtered_main['isoform'].unique())
    logger.info(f"Total transcripts with retained introns: {total_retained_introns}")
    output_file.write(f"Total transcripts with retained introns:\t{total_retained_introns}\n")

    # Initial matches with m6A modifications
    matched_in_m6anet = set(indiv_significant['isoform'].unique()).union(set(indiv_non_significant['isoform'].unique()))
    initial_matches = len(matched_in_m6anet)
    logger.info(f"Initial total retained introns matched in m6anet: {initial_matches}")
    output_file.write(f"Initial total retained introns matched in m6anet:\t{initial_matches}\n")

    # Separate significant and non-significant counts
    passing_threshold = set(indiv_significant['isoform'].unique())
    failing_threshold = set(indiv_non_significant['isoform'].unique()) - passing_threshold  # Ensure no overlap
    
    num_passing = len(passing_threshold)
    num_failing = len(failing_threshold)
    
    logger.info(f"Total retained introns passing modification threshold: {num_passing}")
    logger.info(f"Total retained introns failing modification threshold (subset of matches): {num_failing}")
    output_file.write(f"Total retained introns passing modification threshold:\t{num_passing}\n")
    output_file.write(f"Total retained introns failing modification threshold (subset of matches):\t{num_failing}\n")

    # Calculate the percentages relative to the total retained introns
    perc_matched_in_m6anet = (initial_matches / total_retained_introns) * 100 if total_retained_introns > 0 else 0
    perc_total_passing = (num_passing / total_retained_introns) * 100 if total_retained_introns > 0 else 0
    perc_total_failing = (num_failing / total_retained_introns) * 100 if total_retained_introns > 0 else 0
    
    # Write these new percentage calculations to the output file
    output_file.write(f"Percentage of retained introns identified in m6anet:\t{perc_matched_in_m6anet:.2f}%\n")
    output_file.write(f"Percentage of total retained introns passing threshold:\t{perc_total_passing:.2f}%\n")
    output_file.write(f"Percentage of total retained introns failing threshold:\t{perc_total_failing:.2f}%\n")
    
    logger.info(f"Percentage of retained introns identified in m6anet: {perc_matched_in_m6anet:.2f}%")
    logger.info(f"Percentage of total retained introns passing threshold: {perc_total_passing:.2f}%")
    logger.info(f"Percentage of total retained introns failing threshold: {perc_total_failing:.2f}%")

    # Additional debugging details for further verification
    output_file.write("\n### Debugging Details ###\n")
    output_file.write(f"Retained introns matched: {initial_matches}\n")
    output_file.write(f"Retained introns passing threshold: {num_passing}\n")
    output_file.write(f"Retained introns failing threshold: {num_failing}\n")
    output_file.write(f"Percentage passing threshold (of total retained introns): {perc_total_passing:.2f}%\n")
    output_file.write(f"Percentage failing threshold (of total retained introns): {perc_total_failing:.2f}%\n\n")
    
    logger.info("\n### Debugging Details ###")
    logger.info(f"Retained introns matched: {initial_matches}")
    logger.info(f"Retained introns passing threshold: {num_passing}")
    logger.info(f"Retained introns failing threshold: {num_failing}")
    logger.info(f"Percentage passing threshold (of total retained introns): {perc_total_passing:.2f}%")
    logger.info(f"Percentage failing threshold (of total retained introns): {perc_total_failing:.2f}%\n")

    logger.info("Summary statistics successfully written to the output file.")




def perform_enrichment_analysis(filtered_main, indiv_filtered):
    """
    Perform enrichment analysis to determine if retained introns with m6A modifications
    are overrepresented compared to retained introns without modifications.
    
    Parameters:
    filtered_main (pd.DataFrame): DataFrame with all transcripts with retained introns.
    indiv_filtered (pd.DataFrame): DataFrame with significant modification data from data.indiv_proba.csv.

    Returns:
    tuple: Odds ratio and p-value from the Fisher's Exact Test.
    """
    # Unique counts for retained introns with and without m6A modifications
    total_retained_introns = len(filtered_main['isoform'].unique())
    retained_introns_with_m6a = len(indiv_filtered['isoform'].unique())
    retained_introns_without_m6a = total_retained_introns - retained_introns_with_m6a

    # Construct contingency table for Fisher's Exact Test
    table = [[retained_introns_with_m6a, retained_introns_without_m6a],
             [total_retained_introns - retained_introns_with_m6a, total_retained_introns - retained_introns_without_m6a]]

    # Fisher's Exact Test
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    
    # Output results
    logger.info(f"Enrichment analysis results - Odds Ratio: {odds_ratio}, P-value: {p_value}")
    print(f"Enrichment analysis results - Odds Ratio: {odds_ratio}, P-value: {p_value}")
    return odds_ratio, p_value



def main(args):
    # Load and filter the main data file
    filtered_main = load_and_filter_main_data(args.data)
    mapped_data = map_read_id_to_index(filtered_main, args.index_dir)
    
    # Filter outputs with different matching criteria
    original_sig, original_non_sig = filter_m6anet_files(args.m6anet_dir, mapped_data, match_criteria="original")
    strict_sig, strict_non_sig = filter_m6anet_files(args.m6anet_dir, mapped_data, match_criteria="strict")
    gene_strict_sig, gene_strict_non_sig = filter_m6anet_files(args.m6anet_dir, mapped_data, match_criteria="gene_strict")
    
    # Save all outputs to tab-separated files for verification
    original_sig.to_csv("m6a_retained_intron_data_table_original_significant.tsv", sep='\t', index=False)
    original_non_sig.to_csv("m6a_retained_intron_data_table_original_non_significant.tsv", sep='\t', index=False)
    strict_sig.to_csv("m6a_retained_intron_data_table_strict_significant.tsv", sep='\t', index=False)
    strict_non_sig.to_csv("m6a_retained_intron_data_table_strict_non_significant.tsv", sep='\t', index=False)
    gene_strict_sig.to_csv("m6a_retained_intron_data_table_gene_strict_significant.tsv", sep='\t', index=False)
    gene_strict_non_sig.to_csv("m6a_retained_intron_data_table_gene_strict_non_significant.tsv", sep='\t', index=False)


    # Open a single summary file to log statistics for all filtering levels
    with open("m6a_summary_output.txt", "w") as summary_file:
        
        # Original Output Summary
        summary_file.write("Original Output (gene level match allowed, if trans_ID == novel):\n")
        calculate_summary_stats(filtered_main, original_sig, original_non_sig, summary_file)
        
        # Perform enrichment analysis on significant data only
        odds_ratio, p_value = perform_enrichment_analysis(filtered_main, original_sig)
        summary_file.write(f"Enrichment Analysis - Odds Ratio: {odds_ratio}, P-value: {p_value}\n\n")
        logger.info(f"Original significant enrichment analysis completed - Odds Ratio: {odds_ratio}, P-value: {p_value}")

        # Strict Transcript Match Output Summary
        summary_file.write("########\nSTRICT Transcript to transcript_ID Match Output:\n")
        calculate_summary_stats(filtered_main, strict_sig, strict_non_sig, summary_file)
        
        # Perform enrichment analysis on strict significant data only
        odds_ratio, p_value = perform_enrichment_analysis(filtered_main, strict_sig)
        summary_file.write(f"Enrichment Analysis - Odds Ratio: {odds_ratio}, P-value: {p_value}\n\n")
        logger.info(f"Strict significant enrichment analysis completed - Odds Ratio: {odds_ratio}, P-value: {p_value}")
        
        # Gene-Level Strict Match Output Summary
        summary_file.write("########\nGene-Level Strict Match Output:\n")
        calculate_summary_stats(filtered_main, gene_strict_sig, gene_strict_non_sig, summary_file)
        
        # Perform enrichment analysis on gene strict significant data only
        odds_ratio, p_value = perform_enrichment_analysis(filtered_main, gene_strict_sig)
        summary_file.write(f"Enrichment Analysis - Odds Ratio: {odds_ratio}, P-value: {p_value}\n")
        logger.info(f"Gene-level strict significant enrichment analysis completed - Odds Ratio: {odds_ratio}, P-value: {p_value}")

    # Log summary of row counts for each output for easier comparison
    logger.info(f"Row counts for each output:")
    logger.info(f"Original significant: {original_sig.shape[0]} rows")
    logger.info(f"Original non-significant: {original_non_sig.shape[0]} rows")
    logger.info(f"Strict significant: {strict_sig.shape[0]} rows")
    logger.info(f"Strict non-significant: {strict_non_sig.shape[0]} rows")
    logger.info(f"Gene-level strict significant: {gene_strict_sig.shape[0]} rows")
    logger.info(f"Gene-level strict non-significant: {gene_strict_non_sig.shape[0]} rows")




# Default paths based on your current directory structure
DEFAULT_MAIN_DATA_PATH = "./data/vir_temp_classification.txt"
DEFAULT_INDEX_TO_READS_DIR = "./index_to_reads"
DEFAULT_M6ANET_DIR = "./m6a_net"
DEFAULT_SUMMARY_OUTPUT_PATH = "./m6a_summary_output.txt"

# Setup command-line argument parsing with short options
parser = argparse.ArgumentParser(description="Process m6anet data for retained intron analysis in Arabidopsis.")
parser.add_argument(
    '-d', '--data', 
    type=str, 
    default=DEFAULT_MAIN_DATA_PATH, 
    help="Path to the main data file. Default is './data/vir_temp_classification.txt'."
)
parser.add_argument(
    '-i', '--index_dir', 
    type=str, 
    default=DEFAULT_INDEX_TO_READS_DIR, 
    help="Path to the directory with index to reads files. Default is './index_to_reads'."
)
parser.add_argument(
    '-m', '--m6anet_dir', 
    type=str, 
    default=DEFAULT_M6ANET_DIR, 
    help="Path to the m6anet output directory. Default is './m6a_net'."
)
parser.add_argument(
    '-o', '--output', 
    type=str, 
    default=DEFAULT_SUMMARY_OUTPUT_PATH, 
    help="Path to the summary output file. Default is './m6a_summary_output.txt'."
)



# Uncomment these lines when running as a standalone script
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
