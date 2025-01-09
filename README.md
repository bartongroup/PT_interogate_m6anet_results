# interogate_m6anet_results


Interogate m6Anet results is a Python package for analysing m6A (N6-methyladenosine) modification data and overlaying this with the GTF file. The package provides tools for parsing GTF files, generating transcript coordinates, and querying specific transcript information.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

- Python 3.7 or higher
- Required Python packages: `numpy`, `matplotlib`, `pandas`, `nose2`

### Clone the Repository

```bash
git clone https://github.com/peterthorpe5/interogate_m6anet_results.git
cd interogate_m6anet_results


## Install Required Packages
to test

pip install numpy matplotlib pandas nose nose2

```

# what does this do and how

1) this parses a gft fie (tested) gff3 (not yet tested) and sets up a LARGE dictionary of transcript exon number to
coordinates for the nucleotide sequence.  For exmaple:

AT1G01020.4 exon 2: [283, 284, 285]

In the gtf, the cooridnates are genomic locations, these dont directly help when mapping to the transcriptome. 

```bash
python interogate_m6anet.py

```

### example usage 

```bash

python interogate_m6anet.py 
 --gtf Araprot.11.gtf --m6a */data.site_proba.csv --logfile m6anet_interogator.logfile

```
   The interogate_m6anet.py script processes the initial m6A modification data obtained from m6Anet.

## Additional Processing Scripts


 2) Collect Positions of m6A

After running the `interogate_m6anet.py` script, you need to run `scripts/collect_positions_of_m6a.py`. This script collects and summarizes the positions of m6A modifications identified in the data.

```bash

python scripts/collect_positions_of_m6a.py  --file <infile> --output <name>

--file output/data.site_proba_exon_annotated.tab
```

The input file is the ouput  `data.site_proba_exon_annotated.tab` from  `interogate_m6anet.py` script 

What This Script Does:

-   Extracts positions of m6A modifications from the processed data.
-   Summarizes the positions, providing a detailed list of all identified modification sites.

If you would like to collect all the position from one condtion (i.e. all the rep), then you can use a wild card to collect these

```bash

python scripts/collect_positions_of_m6a.py  --file WT_*/data.site_proba_exon_annotated.tab 
--output WT.results

python scripts/collect_positions_of_m6a.py  --file MUT_*/data.site_proba_exon_annotated.tab 
--output MUT.results

```

3) Compare Positions Between Conditions

With the output from the previous script, run `scripts/compare_positions_between_conditions.py` to compare the positions of m6A modifications between different experimental conditions.

# Compare the results between conditions. 

```bash

python scripts/compare_positions_between_conditions.py --file <infile> --output <name>

```

an example based on the outfiles from above:

```bash

python scripts/compare_positions_between_conditions.py --file WT.results MUT.results 
--output WT_vs_MUT

```

What This Script Does:

- Compares the collected m6A modification positions between different experimental conditions.
- Highlights the differences and similarities in m6A modification sites between the specified conditions.

## test for this script 

```bash

python scripts/compare_positions_between_conditions.py --file tests/test1.sites.summerise 
tests/test2.sites.summerise

```

Note: the WT.results and MUT.results are a collection of all the data from all the reps. 

You can even collect all different permutations of your data as you see fit. Enjoy!! :) 


---

# identify_m6a_in_introns.py

## scripts/identify_m6a_in_introns.py

This is a script to try and identify m6a site that are associated with reads that contain intron as a result of incorrect splicing. 

Script Description: Identifying Retained Introns with m6A Modifications in Arabidopsis - NOTE: Arabidopsis threshold are the defaults. This can be changed if anyone ever reads this. 

`Purpose:`
This script is designed to analyse retained introns within a set of Arabidopsis transcripts to determine if they carry m6A modifications, which are methylations at specific adenosine residues. 

`Overview:`
-d file contains data about isoforms, associated genes, transcripts, and intron subcategories.
Only those transcripts marked with an "intron_retention" subcategory are kept for further analysis.

`Mapping Isoforms to Reads:`
The script then maps each filtered isoform to its associated read_index using files within the index_to_reads directory. This step is essential because read_index is used to match m6A modification data from m6anet outputs.

`Parsing m6A Modification Files:`
Each isoform is analysed to see if it contains m6A modifications based on files in the m6a_net directory, which holds modification probability data.
Each isoform is filtered to identify m6A modifications using two different thresholds: one for individual reads (data.indiv_proba.csv) and one for m6A modification sites (data.site_proba.csv).


`Matching Criteria:`
The script uses three different matching criteria to determine if retained introns with m6A modifications meet specific criteria:
Original: Allows matching at the gene level if the transcript ID is “novel,” meaning an unidentified or newly predicted transcript.
Strict: Requires an exact match between transcript_id and associated_transcript.
Gene Strict: Matches at the gene level without considering the transcript, allowing for more flexibility in finding retained introns with modifications.

`Output and Statistics:`
For each matching criterion, significant and non-significant matches are saved as tab-separated files for later verification and comparison.
Summary statistics for each set are written to an output file. These statistics include the total number of retained introns, those identified with m6A modifications, and those passing or failing the significance threshold.
Percentages are calculated for both significant and non-significant reads, offering insight into the prevalence of m6A modifications within retained introns.

`Enrichment Analysis:`
For each matching criterion, an enrichment analysis is performed using Fisher's Exact Test, comparing the counts of retained introns with and without m6A modifications. This statistical analysis provides an odds ratio and a p-value, helping to determine if m6A modifications are enriched in retained introns.


`Summary of Outputs:`
Filtered Data: Saved as TSV files, differentiating between significant and non-significant m6A modification matches.
Summary Statistics: A  summary output file detailing the counts and percentages of retained introns with m6A modifications.
Enrichment Analysis: Statistical results from Fisher’s Exact Test, offering a p-value and odds ratio for enrichment of m6A modifications in retained introns.


# Documentation for Retained Intron Analysis with m6Anet

This documentation describes the input files, their requirements, and the purpose of the script for analysing retained introns and their m6A methylation modifications using m6Anet outputs.

## Input Files

### 1. **Main Data File** (`-d`)

- **Description**: This file contains classification data for transcripts, including information on intron retention. This is an output file from a flaire/SQANTI3  de novo assembly. Note: A de novo assembly using flaire https://github.com/BrooksLabUCSC/flair is required first (nanopore, and with Illumina is optimal), then SQANTI3 (https://github.com/ConesaLab/SQANTI3) (the output file is actually generated here *_classification.txt). 

- **Format**: Tab-separated values (TSV).
- **Key Columns**:
  - `isoform`: The transcript isoform identifier.
  - `associated_gene`: The gene associated with the transcript.
  - `associated_transcript`: Transcript ID used for further filtering.
  - `subcategory`: Classification of the transcript, such as "intron_retention."
- **Default Path**: `./data/vir_temp_classification.txt`.

---

### 2. **Index to Reads Directory** (`-i`)

- **Description**: Directory containing mappings of transcript isoform IDs to `read_index` values, which are required for aligning the main data file to m6Anet outputs. This file is generated using f5c (https://hasindu2008.github.io/f5c/docs/commands) with the option --summary FILE. F5c is used to generate the event align file for m6anet (and is faster)
- **File Requirements**:
  - **Format**: Tab-separated values (TSV).
  - **Key Columns**:
    - `read_name`: The isoform identifier.
    - `read_index`: The corresponding index required for merging with m6Anet data.
  - File names must end with `.tsv`.
- **Default Path**: `./index_to_reads`.  can have a directory of files, it will use them all. 

---

### 3. **m6Anet Output Directory** (`-m`)

- **Description**: Directory containing output files generated by m6Anet, which include probability scores for RNA modification. This is the standard ouput. 
- **Structure**:
  - Each subdirectory corresponds to a processed sample and contains two key files:
    1. `data.indiv_proba.csv`: Contains individual modification probabilities.
    2. `data.site_proba.csv`: Contains site-level probabilities for modifications.
- **Key Columns in `data.indiv_proba.csv`**:
  - `read_index`: Index to align modifications with transcript data.
  - `transcript_id`: Transcript ID for strict matching.
  - `probability_modified`: Probability score for modification.
- **Default Path**: `./m6a_net`.

---

## Input Validation

1. **Main Data File**:
   - Ensure the `subcategory` column contains values such as "intron_retention" for filtering.
   - Must be a valid TSV file with the specified columns.

2. **Index to Reads Directory**:
   - All files must be non-empty and contain the required columns (`read_name`, `read_index`).
   - Skips empty files and logs warnings for missing data.

3. **m6Anet Directory**:
   - Each subdirectory must include `data.indiv_proba.csv`.
   - Logs warnings if the required columns (`read_index`, `probability_modified`) are missing.

---

## Outputs

### 1. Filtered Data Files

The retained intron usually cannot be assigned speicifcally to a transcript, as these things are usually messed up and with transcript does it most closely match???!!! So, really this can only be confidently reported on a per gene basis, thus the "gene_strict_significant" ouput file is the most useful for you. 

- **Significant Data**: Contains retained introns with m6A modifications passing the probability threshold.
- **Non-Significant Data**: Contains retained introns with m6A modifications failing the probability threshold.
- **Files Generated**:
  - `m6a_retained_intron_data_table_original_significant.tsv`
  - `m6a_retained_intron_data_table_original_non_significant.tsv`
  - `m6a_retained_intron_data_table_strict_significant.tsv`
  - `m6a_retained_intron_data_table_strict_non_significant.tsv`
  - `m6a_retained_intron_data_table_gene_strict_significant.tsv`
  - `m6a_retained_intron_data_table_gene_strict_non_significant.tsv`

  # Explanation of Matching Criteria: Original, Strict, and Gene-Strict

The terms **strict** and **original** refer to different matching criteria used in the script to align the data between the main dataset and the m6Anet outputs. Here's a detailed explanation of what each term means:

---

## **Original Matching**
- **Definition**: The least restrictive matching criterion, allowing flexibility in aligning datasets.
- **Key Features**:
  - Matches are primarily based on the `read_index` column.
  - It allows for matches between `associated_transcript` in the main data and `transcript_id` in the m6Anet output even if:
    - The `associated_transcript` is labelled as `"novel"`, and the `associated_gene` matches the gene part of the `transcript_id` (before the dot).
    - This accommodates novel transcripts or situations where the transcript ID may not be a direct match but can be inferred from the gene.
- **Purpose**: To maximise the number of matches, even if they are less stringent. This is useful when novel transcripts or partial data are present.

---

## **Strict Matching**
- **Definition**: A more restrictive criterion that ensures matches are highly specific and accurate.
- **Key Features**:
  - Matches are based on `read_index` and `associated_transcript` aligning exactly with `transcript_id` in the m6Anet output.
  - No allowance is made for `"novel"` transcript IDs or partial matches.
  - Ensures only well-defined and explicitly matching transcripts are included.
- **Purpose**: To provide high-confidence matches, minimising the inclusion of potentially ambiguous or incorrectly assigned transcripts.

---

## **Gene-Strict Matching**
- **Definition**: A criterion that focuses on matching at the gene level rather than the transcript level.
- **Key Features**:
  - Matches are based on `read_index` and ensure that `associated_gene` aligns with the gene part of `transcript_id` (before the dot).
  - Less restrictive than strict matching at the transcript level but ensures consistency at the gene level.
- **Purpose**: To account for cases where transcript-level matching may be ambiguous, but gene-level relationships are still valid.

---

## **Summary of Differences**
| **Criterion**      | **Matching Basis**                       | **Flexibility**                                                                                     | **Purpose**                                              |
|---------------------|------------------------------------------|-----------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| **Original**        | `read_index` and `associated_transcript` | Allows `"novel"` transcripts and partial matches (e.g., gene-level matches for novel transcripts).  | Maximises matches, useful for partial or novel data.     |
| **Strict**          | `read_index` and `associated_transcript` | Requires exact matches between `associated_transcript` and `transcript_id`.                        | Provides high-confidence matches.                       |
| **Gene-Strict**     | `read_index` and `associated_gene`       | Matches `associated_gene` with the gene part of `transcript_id` (before the dot).                  | Focuses on gene-level consistency, allows flexibility.  |

This flexibility is provided to adapt the script to different levels of stringency depending on the confidence or quality of the data being analysed.


### 2. Summary Output File

- **Description**: A text file summarising the analysis, including:
  - Total transcripts with retained introns.
  - Matches in m6Anet significant and non-significant categories.
  - Percentages of retained introns with m6A modifications.
  - Results of enrichment analysis.
- **Generated File**: `m6a_summary_output.txt`.

---

## Additional Notes

### Matching Criteria

1. **Original**:
   - Matches based on `read_index` and allows for "novel" transcripts to match genes.
2. **Strict**:
   - Matches strictly based on `read_index` and `associated_transcript` aligning to `transcript_id`.
3. **Gene-Strict**:
   - Matches based on `read_index` and verifies that `associated_gene` aligns to the gene part of `transcript_id`.

### Logging

- Logs are written to both the console and a log file (`intron_retention_m6anet_methylated.log`).
- Provides detailed information on filtering, merging, and statistical analysis steps.

---

## Default Paths

- **Main Data File**: `./data/vir_temp_classification.txt` -  that was my testing file. 
- **Index Directory**: `./index_to_reads`
- **m6Anet Directory**: `./m6a_net`
- **Summary Output File**: `./m6a_summary_output.txt`
