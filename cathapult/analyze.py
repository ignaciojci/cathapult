import pandas as pd
import re
from pathlib import Path

def load_reference_data(data_dir):
    # Load CATH names
    cath_names = []
    with open(data_dir / "cath-names.txt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                left, name = line.split(":", 1)
                cath_id = left.strip().split()[-2]
                cath_names.append((cath_id, name.strip()))
    df_cath_names = pd.DataFrame(cath_names, columns=["CATH_ID", "NAME"])

    # Load CATH superfamily list (tab-delimited)
    df_cath_super = pd.read_csv(data_dir / "cath-superfamily-list.txt", sep="\t", dtype=str)

    return df_cath_names, df_cath_super

import pandas as pd

def extract_domain_levels(df):
    # Check column exists
    if "domain" not in df.columns:
        raise ValueError("Column 'domain' not found in DataFrame")

    # Work on a copy to avoid modifying original unintentionally
    df = df.copy()

    # Ensure string type and handle None/NaN safely
    domain_series = df["domain"].fillna("").astype(str)

    # Split safely
    domain_split = domain_series.str.split(".", expand=True)

    # First level (always safe if string is not empty)
    df["domain.first.level"] = domain_split[0].replace("", pd.NA)

    # Two levels (only if second level exists)
    if domain_split.shape[1] >= 2:
        df["domain.two.levels"] = (
            domain_split[[0, 1]]
            .apply(lambda x: ".".join(filter(None, x)), axis=1)
            .replace("", pd.NA)
        )
    else:
        df["domain.two.levels"] = pd.NA

    # Three levels (only if third level exists)
    if domain_split.shape[1] >= 3:
        df["domain.three.levels"] = (
            domain_split[[0, 1, 2]]
            .apply(lambda x: ".".join(filter(None, x)), axis=1)
            .replace("", pd.NA)
        )
    else:
        df["domain.three.levels"] = pd.NA

    return df

def count_domains(df, col, unique_gene=False):
    sub = df[["gene", col]].dropna()
    if unique_gene:
        sub = sub.drop_duplicates()
    counts = sub[col].value_counts().reset_index()
    counts.columns = ["domain", "count"]
    return counts

def annotate_domains(counts_df):
    # Load reference data
    data_dir = Path(__file__).parent / "data"
    cath_names, cath_super = load_reference_data(data_dir)
    
    # Deduplicate CATH names to avoid InvalidIndexError
    cath_names_unique = cath_names.drop_duplicates(subset="CATH_ID", keep="first")

    # First: match using cath-names.txt
    annotated = counts_df.copy()
    annotated["domain.name"] = annotated["domain"].map(
        cath_names_unique.set_index("CATH_ID")["NAME"]
    )

    # Second: fill remaining names using cath-superfamily-list.txt
    annotated["domain.name"] = annotated["domain.name"].fillna(
        cath_super.set_index("# CATH_ID")["NAME"]
    )

    return annotated

def deep_annotate_domains(df):
    # Load reference data
    data_dir = Path(__file__).parent / "data"
    cath_names, cath_super = load_reference_data(data_dir)
    
    # Deduplicate CATH names to avoid InvalidIndexError
    cath_names_unique = cath_names.drop_duplicates(subset="CATH_ID", keep="first")

    df = extract_domain_levels(df)
    
    # First: match using cath-names.txt
    annotated = df.copy()
    annotated["domain.name"] = annotated["domain"].map(
        cath_names_unique.set_index("CATH_ID")["NAME"]
    )

    # Second: fill remaining names using cath-superfamily-list.txt
    annotated["domain.name"] = annotated["domain.name"].fillna(
        cath_super.set_index("# CATH_ID")["NAME"]
    )
    
    annotated["domain.first.level.name"] = annotated["domain.first.level"].map(
        cath_names_unique.set_index("CATH_ID")["NAME"]
    )

    annotated["domain.second.level.name"] = annotated["domain.two.levels"].map(
        cath_names_unique.set_index("CATH_ID")["NAME"]
    )

    annotated["domain.third.level.name"] = annotated["domain.three.levels"].map(
        cath_names_unique.set_index("CATH_ID")["NAME"]
    )
    
    return annotated

def analyze_ted_summary(input_tsv, output_tsv=None):
    df = pd.read_csv(input_tsv, sep="\t")
    
    df = df[df["cath_label"] != "-"]
    df["gene"] = df["ted_id"].str.extract(r"AF-(.*)-F1")
    df["domain"] = df["cath_label"]

    df = extract_domain_levels(df)

    domain_cols = ["domain", "domain.first.level", "domain.two.levels", "domain.three.levels"]
    results = []

    for col in domain_cols:
        for dedup in [False, True]:
            counts = count_domains(df, col, unique_gene=dedup)
            counts["domain.type"] = col
            counts["deduplicated"] = "deduped" if dedup else ""
            results.append(counts)

    combined = pd.concat(results, ignore_index=True)

    annotated = annotate_domains(combined)
    
    if output_tsv:
        annotated.to_csv(output_tsv, sep="\t", index=False)
        print(f"Saved annotated domain counts to {output_tsv}")
    return annotated