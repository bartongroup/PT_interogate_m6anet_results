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

python scripts/collect_positions_of_m6a.py  --file WT_*/data.site_proba_exon_annotated.tab --output WT.results

python scripts/collect_positions_of_m6a.py  --file MUT_*/data.site_proba_exon_annotated.tab --output MUT.results
```

3) Compare Positions Between Conditions

With the output from the previous script, run `scripts/compare_positions_between_conditions.py` to compare the positions of m6A modifications between different experimental conditions.

# Compare the results between conditions. 

```bash

python scripts/compare_positions_between_conditions.py --file <infile> --output <name>

```

an example based on the outfiles from above:

```bash

python scripts/compare_positions_between_conditions.py --file WT.results MUT.results --output WT_vs_MUT

```

What This Script Does:

- Compares the collected m6A modification positions between different experimental conditions.
- Highlights the differences and similarities in m6A modification sites between the specified conditions.

Note: the WT.results and MUT.results are a collection of all the data from all the reps. 

You can even collect all different permutations of your data as you see fit. Enjoy!! :) 

