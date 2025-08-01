# cathapult

**cathapult** is a lightweight CLI tool and Python package for pulling, filtering, analyzing, and comparing CATH-TED domain annotations from UniProt IDs.
It streamlines the extraction of domain summaries, counts domain occurrences, and performs enrichment analysis between two groups.

---

## Features

* **Fetch** domain summaries for specific UniProt IDs from CATH-TED.
* **Filter** a large CATH-TED domain summary file (`.tsv` or `.tsv.gz`) by a set of UniProt IDs and an optional keyword.
* **Count** CATH domains at multiple hierarchy levels.
* **Annotate** domain IDs with human-readable names.
* **Compare** feature enrichment between two groups using odds ratio analysis and visualization. 📊
* Works as both a **command-line tool** and a **Python API**.

---

## Installation

```bash
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install .
```

---

## Quick Start

### 1. Fetch domain summaries from UniProt IDs

Create a text file (`target_uniprot_ids.txt`) with one UniProt ID per line, then:

```bash
cathapult fetch target_uniprot_ids.txt output/summary.tsv
```

*The output directory is created automatically if it doesn’t exist.*

---

### 2. Filter a large CATH-TED domain summary file

If you already have a large domain summary file (e.g., from bulk download),
you can filter it by your target UniProt IDs and an optional keyword (default: `"sapiens"`).

```bash
cathapult filter target_uniprot_ids.txt --domain_summary_file path/to/domain_summary.tsv.gz --keyword sapiens --output_file filtered.tsv
```

💡 **Tip:** You can skip `--domain_summary_file` if you set the environment variable:

```bash
export DOMAIN_SUMMARY_FILE=/path/to/domain_summary.tsv.gz   # macOS/Linux
setx DOMAIN_SUMMARY_FILE "C:\path\to\domain_summary.tsv.gz" # Windows (new shells)
```

Then simply run:

```bash
cathapult filter target_uniprot_ids.txt --keyword sapiens --output_file filtered.tsv
```

---

### 3. Analyze a filtered summary

Count and annotate the domains from a fetched or filtered summary file:

```bash
cathapult analyze filtered.tsv output/domain_counts.tsv
```

---

### 4. Compare feature enrichment between two groups

Use the `odds-ratio` command to identify enriched domains between two datasets:

```bash
cathapult odds-ratio group1.tsv group2.tsv \
    --output_plot enrichment_plot.png \
    --output_table enrichment_table.tsv
```

Options:

* `--unique_features` → Count duplicated features within each protein only once.
* `--filter_column` and `--filter_values` → Restrict analysis to specific categories.
* `--alpha` → Significance threshold for coloring plot points.

---

## Python API

You can also run the same operations programmatically:

```python
from cathapult.fetcher import fetch_ted_summary
from cathapult.fetcher import filter_ted_summary
from cathapult.analyze import analyze_ted_summary

# Fetch summaries
ids = ["P12345", "Q67890"]
all_data = [row for uid in ids for row in fetch_ted_summary(uid)]
  
# Filter from a bulk domain summary file
df_filtered = filter_ted_summary(
    target_ids=set(ids),
    ted_db_gz="domain_summary.tsv.gz",
    filter_keyword="sapiens"
)

# Analyze
df_counts = analyze_ted_summary("filtered.tsv")
```

---

## Development

Clone and install in editable mode:

```bash
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install -e .
```

---

## Acknowledgments

* [CATH Database](https://www.cathdb.info/)
* TED Databases
* UniProt
