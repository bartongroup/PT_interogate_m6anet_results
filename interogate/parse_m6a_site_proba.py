#!/usr/bin/env python3
#
# interogate_m6anet.py 

import pandas as pd


def identify_methylated_sites(m6a_site_proba, threshold=0.9):
    """
    Identify methylated sites with probability greater than the threshold.

    Parameters:
    m6a_site_proba (str): Path to the CSV file.
    threshold (float): Probability threshold to consider for methylation prediction.

    Returns:
    pd.DataFrame: DataFrame containing transcript ID and positions of methylated sites above the threshold.
    """
    # Load the CSV file into a DataFrame
    df = pd.read_csv(m6a_site_proba)

    # Ensure the necessary columns exist
    required_columns = ['transcript_id', 'transcript_position', 'probability_modified']
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"The input file must contain the following columns: {required_columns}")

    # Filter rows where the probability is greater than the threshold
    methylated_sites = df[df['probability_modified'] > threshold]

    # Select and return the relevant columns
    result = methylated_sites[['transcript_id', 'transcript_position']]
    return result





