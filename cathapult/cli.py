import argparse
import os
import pandas as pd
from .fetcher import fetch_ted_summary
from .analyze import analyze_ted_summary

def cli_fetch(args):
    # Determine output path
    if args.output_file:
        output_file = args.output_file
    else:
        base = os.path.splitext(os.path.basename(args.input_file))[0]
        output_file = f"{base}.tsv"

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Read input UniProt IDs
    with open(args.input_file) as f:
        uniprot_ids = [line.strip() for line in f if line.strip()]

    print(f"Fetching domain summaries for {len(uniprot_ids)} UniProt IDs...")

    all_data = []
    for uid in uniprot_ids:
        data = fetch_ted_summary(uid, delay=0.1)
        all_data.extend(data)

    if all_data:
        df = pd.DataFrame(all_data)
        df.to_csv(output_file, sep="\t", index=False)
        print(f"Results saved to: {output_file}")
    else:
        print("No data fetched.")

def cli_analyze(args):
    analyze_ted_summary(args.input_file, args.output_file)

def main():
    parser = argparse.ArgumentParser(description="CATH-TED domain fetcher and analyzer")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Fetch subcommand
    fetch_parser = subparsers.add_parser("fetch", help="Fetch CATH-TED domain summaries from UniProt IDs")
    fetch_parser.add_argument("input_file", help="Text file with UniProt IDs (one per line)")
    fetch_parser.add_argument("output_file", nargs="?", help="Optional output TSV path")
    fetch_parser.set_defaults(func=cli_fetch)

    # Analyze subcommand
    analyze_parser = subparsers.add_parser("analyze", help="Count and annotate domain types from a fetched TSV")
    analyze_parser.add_argument("input_file", help="TSV file from `fetch` command")
    analyze_parser.add_argument("output_file", help="Output TSV file with domain counts and annotations")
    analyze_parser.set_defaults(func=cli_analyze)

    # Parse and dispatch
    args = parser.parse_args()
    args.func(args)