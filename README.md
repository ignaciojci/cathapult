# cathapult

**cathapult** is a lightweight CLI tool and Python package for pulling and analyzing CATH-TED domain annotations from Uniprot IDs. It streamlines the extraction of UniProt-specific domain summaries, counts domain occurrences, and annotates them with human-readable names from official CATH resources.

---

## Features

- Filter domain summaries for specific UniProt IDs.
- Count CATH domains at multiple hierarchy levels.
- Annotate domain IDs with human-readable names.
- Command-line and Python API support.

---

## Installation

```
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install .
```


⸻

Quick Start

1. Filter domain summaries by UniProt ID

```
cathapult fetch \
  --input target_uniprot_ids.txt \
  --summary ted_365m.domain_summary.cath.globularity.taxid.tsv.gz \
  --output output/filtered_summary.tsv
```

Output directory is created automatically if it doesn’t exist.

⸻

2. Analyze filtered summary

```
cathapult analyze \
  --input output/filtered_summary.tsv \
  --output output/count_cath_domains.tsv
```

This will:
	•	Count domains at full and truncated CATH ID levels.
	•	Optionally deduplicate by gene.
	•	Annotate using cath-names.txt and cath-superfamily-list.txt.

⸻

Input File Format

UniProt ID list

Plain text file with one UniProt ID per line:
```
Q9Y2T3
P38398
```

TED domain summary

Expected to be a *.tsv file as output from CATH-TED. Can be gzipped.

⸻

Python API

You can also use cathapult programmatically:

from cathapult.analyze import analyze_ted_summary

df = analyze_ted_summary("filtered_summary.tsv")
df.to_csv("count_cath_domains.tsv", sep="\t", index=False)

⸻

Development

Clone and install in editable mode:

```
git clone https://github.com/ignaciojci/cathapult.git
cd cathapult
pip install -e .
```
⸻

Acknowledgments
- CATH and TED Database
- UniProt