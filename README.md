# cathapult

**cathapult** is a lightweight CLI tool and Python package for pulling, filtering, analyzing, and comparing CATH-TED domain annotations from UniProt IDs.

It provides two main workflows:

1.  **Direct fetching** from APIs for small sets of proteins.
2.  **High-performance local querying** against a pre-built DuckDB database for handling massive, gigabyte-scale domain summary files with response times in a few seconds.

-----

## Features

  * **Fetch** domain summaries for specific UniProt IDs from CATH-TED.
  * **Build** a local, compressed DuckDB database from a massive domain summary file.
  * **Query** the local database to filter millions of rows by UniProt IDs and keywords in a few seconds.
  * **Filter** large summary files directly (in-memory) without a database.
  * **Count** CATH domains at multiple hierarchy levels.
  * **Annotate** domain IDs with human-readable names.
  * **Compare** feature enrichment between two groups using odds ratio analysis and visualization. 📊
  * Works as both a **command-line tool** and a **Python API**.

-----

## Installation

First, install the `duckdb` dependency. Using a binary is recommended for speed and compatibility.

```bash
# It is recommended to install the duckdb binary first
pip install --only-binary :all: duckdb
```

Then, install **cathapult** from the git repository:

```bash
pip install git+https://github.com/ignaciojci/cathapult.git
```

-----

## High-Performance Workflow (Recommended for Large Files)

This workflow is designed for speed when you are working with a large domain summary file, like the 365 million row file from Zenodo. You first build a database once, and then you can query it almost instantly many times.

### 1\. Download the Full Domain Summary

Download a comprehensive domain summary file, such as [`ted_365m.domain_summary.cath.globularity.taxid.tsv.gz`](https://www.google.com/search?q=%5Bhttps://zenodo.org/records/13908086%5D\(https://zenodo.org/records/13908086\)).

### 2\. Build the Local Database

Convert the gzipped TSV into a highly efficient DuckDB database. This is a one-time setup step that makes all future filtering operations extremely fast.

```bash
# This will create ted_365m.domain_summary.cath.globularity.taxid.tsv.duckdb
cathapult setup-db path/to/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz
```

You can use `--overwrite` to rebuild an existing database.

### 3\. Query the Database

Now, use the `query` command to filter the database by your UniProt IDs. This will be significantly faster than reading the gzipped text file each time.

Create a text file (`target_uniprot_ids.txt`) with one UniProt ID per line, then:

```bash
cathapult query path/to/your.duckdb target_uniprot_ids.txt --output_file filtered_results.tsv
```

You can also add a keyword filter:

```bash
cathapult query path/to/your.duckdb target_uniprot_ids.txt --keyword Human --output_file human_filtered.tsv
```

-----

## Alternative Workflows

### Fetching Directly from UniProt

If you only have a few UniProt IDs and don't need the full database, you can fetch their summaries directly.

```bash
cathapult fetch target_uniprot_ids.txt output/summary.tsv
```

*The output directory is created automatically if it doesn’t exist.*

### Filtering Without a Database

If you want to filter the large `.tsv.gz` file without the initial `setup-db` step, you can use the `filter` command. This reads the file in-memory and can be slower for repeated operations.

```bash
cathapult filter target_uniprot_ids.txt path/to/domain_summary.tsv.gz --keyword Human --output_file filtered.tsv
```

💡 **Tip:** The `domain_summary.tsv.gz` path is a positional argument. You can also skip providing it in the command if you set an environment variable:

```bash
# macOS/Linux
export DOMAIN_SUMMARY_FILE="/path/to/domain_summary.tsv.gz"

# Windows (new shells)
setx DOMAIN_SUMMARY_FILE "C:\path\to\domain_summary.tsv.gz" 
```

Then simply run:

```bash
cathapult filter target_uniprot_ids.txt --keyword Human --output_file filtered.tsv
```

-----

## Downstream Analysis

Once you have a filtered TSV file from `query`, `fetch`, or `filter`, you can perform downstream analyses.

### Analyze a Filtered Summary

Count and annotate the domains from a summary file:

```bash
cathapult analyze filtered_results.tsv output/domain_counts.tsv
```

### Compare Feature Enrichment Between Two Groups

Use the `odds-ratio` command to identify enriched domains between two datasets:

```bash
cathapult odds-ratio group1.tsv group2.tsv \
    --output_plot enrichment_plot.png \
    --output_table enrichment_table.tsv
```

**Options:**

  * `--unique_features`: Count duplicated features within each protein only once.
  * `--filter_column` and `--filter_values`: Restrict analysis to specific categories.
  * `--alpha`: Significance threshold for coloring plot points.

-----

## Python API

You can also run the same operations programmatically:

```python
from cathapult.fetcher import fetch_ted_summary, filter_ted_summary
from cathapult.analyze import analyze_ted_summary
from cathapult.db import create_db, query_by_uniprot_ids

# --- Option 1: Fetching and Filtering ---
ids = ["P12345", "Q67890"]

# Fetch summaries directly
all_data = [row for uid in ids for row in fetch_ted_summary(uid)]
  
# Filter from a bulk domain summary file (in-memory)
df_filtered = filter_ted_summary(
    target_ids=set(ids),
    ted_db_gz="domain_summary.tsv.gz",
    filter_keyword="Human"
)

# --- Option 2: Using the Database API ---

# Create the database once
db_path = create_db(tsv_gz_path="domain_summary.tsv.gz")

# Query the database
df_query_results = query_by_uniprot_ids(
    uniprot_ids=ids,
    db_path=db_path,
    keyword="Human"
)

# --- Downstream Analysis ---
df_counts = analyze_ted_summary("filtered_results.tsv")
```

-----

## Development

Clone and install in editable mode:

```bash
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install -e .
```

-----

## Acknowledgments

  * [CATH Database](https://www.cathdb.info/)
  * TED Databases
  * UniProt