import argparse
import os
import pandas as pd
from .fetcher import fetch_ted_summary
from .analyze import analyze_ted_summary
from .enrichment import calculate_odds_ratio, plot_odds_ratio

def cli_fetch(args):
    """Handler for the 'fetch' command."""
    if args.output_file:
        output_file = args.output_file
    else:
        base = os.path.splitext(os.path.basename(args.input_file))[0]
        output_file = f"{base}.tsv"

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    with open(args.input_file) as f:
        uniprot_ids = [line.strip() for line in f if line.strip()]

    print(f"Fetching domain summaries for {len(uniprot_ids)} UniProt IDs...")
    all_data = [data for uid in uniprot_ids for data in fetch_ted_summary(uid, delay=0.1)]

    if all_data:
        pd.DataFrame(all_data).to_csv(output_file, sep="\t", index=False)
        print(f"Results saved to: {output_file}")
    else:
        print("No data fetched.")

def cli_analyze(args):
    """Handler for the 'analyze' command."""
    analyze_ted_summary(args.input_file, args.output_file)

def cli_odds_ratio(args):
    """Handler for the 'odds-ratio' command."""
    print(f"Loading group 1 data from: {args.group1_file}")
    df1 = pd.read_csv(args.group1_file, sep='\t')
    
    print(f"Loading group 2 data from: {args.group2_file}")
    df2 = pd.read_csv(args.group2_file, sep='\t')

    if args.filter_column and args.filter_values:
        print(f"Filtering data on column '{args.filter_column}' with values: {args.filter_values}")
        df1 = df1[df1[args.filter_column].isin(args.filter_values)]
        df2 = df2[df2[args.filter_column].isin(args.filter_values)]
    
    print(f"Calculating odds ratio...")
    results = calculate_odds_ratio(df1, df2, args.unique_features)

    if not results.empty:
        if args.output_table:
            results.to_csv(args.output_table, sep="\t", index=False)
        plot_odds_ratio(results, args.output_plot, args.alpha)
    else:
        print("No features found to plot.")

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(description="CATH-TED domain fetcher and analyzer")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # Fetch subcommand
    fetch_parser = subparsers.add_parser("fetch", help="Fetch CATH-TED domain summaries from UniProt IDs")
    fetch_parser.add_argument("input_file", help="Text file with UniProt IDs (one per line)")
    fetch_parser.add_argument("output_file", nargs="?", help="Optional output TSV path (default: <input_name>.tsv)")
    fetch_parser.set_defaults(func=cli_fetch)

    # Analyze subcommand
    analyze_parser = subparsers.add_parser("analyze", help="Count and annotate domain types from a fetched TSV")
    analyze_parser.add_argument("input_file", help="TSV file from `fetch` command")
    analyze_parser.add_argument("output_file", help="Output TSV file with domain counts and annotations")
    analyze_parser.set_defaults(func=cli_analyze)

    # NEW: Odds-ratio subcommand
    or_parser = subparsers.add_parser(
        "odds-ratio",
        help="Compute and visualize odds ratio of enriched features between two groups."
    )
    or_parser.add_argument("group1_file", help="Path to the TSV file for group 1.")
    or_parser.add_argument("group2_file", help="Path to the TSV file for group 2.")
    or_parser.add_argument(
        "--output_plot", default="odds_ratio_plot.png",
        help="Name of the output plot file (default: 'odds_ratio_plot.png')."
    )
    or_parser.add_argument(
        "--output_table", default="odds_ratio_table.tsv",
        help="Name of the output table file (default: 'odds_ratio_table.tsv')."
    )
    or_parser.add_argument(
        "--unique_features", action='store_true',
        help="Count duplicated features within each protein only once."
    )
    or_parser.add_argument("--filter_column", help="Column to filter data by (optional).")
    or_parser.add_argument(
        "--filter_values", nargs='+',
        help="Values to include for the filter column (requires --filter_column)."
    )
    or_parser.add_argument(
        "--alpha", type=float, default=0.05,
        help="Significance level (alpha) for coloring plot points (default: 0.05)."
    )
    or_parser.set_defaults(func=cli_odds_ratio)

    # Parse and dispatch
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()