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



# in script: identify_m6a_in_introns.py

The input is the ouput from Flaire/ sqantii ... 


Script Description: Identifying Retained Introns with m6A Modifications in Arabidopsis
Purpose:
This script is designed to analyse retained introns within a set of Arabidopsis transcripts to determine if they carry m6A modifications, which are methylations at specific adenosine residues. Understanding the prevalence and distribution of these m6A modifications within retained introns can provide insight into gene regulation, particularly in the context of RNA processing and stability. By categorising these modifications based on their significance and matching criteria, this script enables researchers to explore m6A enrichment and compare between sets of introns with and without modifications.

Overview:
Input Parsing and Filtering:
The script reads an initial classification file (vir_temp_classification.txt by default), which contains data about isoforms, associated genes, transcripts, and intron subcategories.
Only those transcripts marked with an "intron_retention" subcategory are kept for further analysis.
Mapping Isoforms to Reads:
The script then maps each filtered isoform to its associated read_index using files within the index_to_reads directory. This step is essential because read_index is used to match m6A modification data from m6anet outputs.
Parsing m6A Modification Files:
Each isoform is analysed to see if it contains m6A modifications based on files in the m6a_net directory, which holds modification probability data.
Each isoform is filtered to identify m6A modifications using two different thresholds: one for individual reads (data.indiv_proba.csv) and one for m6A modification sites (data.site_proba.csv).
Matching Criteria:
The script uses three different matching criteria to determine if retained introns with m6A modifications meet specific criteria:
Original: Allows matching at the gene level if the transcript ID is “novel,” meaning an unidentified or newly predicted transcript.
Strict: Requires an exact match between transcript_id and associated_transcript.
Gene Strict: Matches at the gene level without considering the transcript, allowing for more flexibility in finding retained introns with modifications.
Output and Statistics:
For each matching criterion, significant and non-significant matches are saved as tab-separated files for later verification and comparison.
Detailed summary statistics for each set are written to an output file. These statistics include the total number of retained introns, those identified with m6A modifications, and those passing or failing the significance threshold.
Percentages are calculated for both significant and non-significant reads, offering insight into the prevalence of m6A modifications within retained introns.
Enrichment Analysis:
For each matching criterion, an enrichment analysis is performed using Fisher's Exact Test, comparing the counts of retained introns with and without m6A modifications. This statistical analysis provides an odds ratio and a p-value, helping to determine if m6A modifications are enriched in retained introns.
Logging and Debugging:
The script incorporates detailed logging for easy tracking of each step, enabling users to follow the parsing, filtering, and merging processes.
Debugging statements output additional details on the number of rows processed, retained introns with and without modifications, and percentage calculations. This is critical for verifying the accuracy of calculations and understanding the dataset structure.
Why This Script is Useful:
This script is specifically designed to support researchers studying the role of RNA modifications in plant gene regulation. m6A modifications have been linked to numerous biological processes, including alternative splicing, mRNA decay, and translation efficiency. By focusing on retained introns, the script provides a targeted approach to explore non-canonical splicing events and the potential regulatory impact of m6A modifications on these sequences. Moreover, the flexibility in matching criteria allows users to tailor the analysis based on gene models, novel transcripts, or strict transcript annotations, making it adaptable to a range of research questions in RNA biology and plant genomics.

Summary of Outputs:
Filtered Data: Saved as TSV files, differentiating between significant and non-significant m6A modification matches.
Summary Statistics: A comprehensive summary output file detailing the counts and percentages of retained introns with m6A modifications.
Enrichment Analysis: Statistical results from Fisher’s Exact Test, offering a p-value and odds ratio for enrichment of m6A modifications in retained introns.
Detailed Logs: Logs capture each step of the process, highlighting successes and any issues encountered with the data files.
This script provides a robust approach to assess RNA methylation within specific intron retention events, helping researchers to understand the modification landscape and its possible functional implications in Arabidopsis and beyond.

