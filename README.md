# cathapult

**cathapult** is a lightweight CLI tool and Python package for pulling, analyzing, and comparing CATH-TED domain annotations from Uniprot IDs. It streamlines the extraction of domain summaries, counts domain occurrences, and performs enrichment analysis between two groups.

-----

## Features

  - Fetch domain summaries for specific UniProt IDs from CATH-TED.
  - Count CATH domains at multiple hierarchy levels.
  - Annotate domain IDs with human-readable names.
  - Compare feature enrichment between two groups using odds ratio analysis and visualization. 📊
  - Command-line and Python API support.

-----

## Installation

```bash
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install .
```

-----

## Quick Start

### 1\. Fetch domain summaries by UniProt ID

First, create a text file (`target_uniprot_ids.txt`) with one UniProt ID per line. Then, run `fetch` to get the domain data.

```bash
cathapult fetch target_uniprot_ids.txt output/filtered_summary.tsv
```

*The output directory is created automatically if it doesn’t exist.*

-----

### 2\. Analyze a filtered summary

Next, you can count and annotate the domains from a fetched summary file.

```bash
cathapult analyze output/filtered_summary.tsv output/count_cath_domains.tsv
```

*This counts domains and annotates them with names.*

-----

### 3\. Compare feature enrichment between two groups

Finally, use the `odds-ratio` command to compare two summary files (e.g., case vs. control) and identify enriched domains.

```bash
cathapult odds-ratio group1_summary.tsv group2_summary.tsv \
    --output_plot enrichment_plot.png \
    --output_table enrichment_table.tsv
```

This command calculates the odds ratio for each domain, generates a plot visualizing the enrichment, and saves the full statistical results to a TSV file.

-----

## Python API

You can also use `cathapult` programmatically:

```python
from cathapult.analyze import analyze_ted_summary

# Analyze a single summary file
df = analyze_ted_summary("filtered_summary.tsv")
df.to_csv("count_cath_domains.tsv", sep="\t", index=False)
```

-----

## Development

Clone the repository and install it in editable mode to start developing.

```bash
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install -e .
```

-----

## Acknowledgments

  - CATH and TED Databases
  - UniProt