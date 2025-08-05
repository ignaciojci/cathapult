import duckdb
import pandas as pd
from pathlib import Path
import gzip
import re
import os
from typing import List, Optional

DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS domain_summary (
    ted_id TEXT,
    md5_domain TEXT,
    consensus_level TEXT,
    chopping TEXT,
    nres_domain INTEGER,
    num_segments INTEGER,
    plddt REAL,
    num_helix_strand_turn INTEGER,
    num_helix INTEGER,
    num_strand INTEGER,
    num_helix_strand INTEGER,
    num_turn INTEGER,
    proteome_id TEXT,
    cath_label TEXT,
    cath_assignment_level TEXT,
    cath_assignment_method TEXT,
    packing_density REAL,
    norm_rg REAL,
    tax_common_name TEXT,
    tax_scientific_name TEXT,
    tax_lineage TEXT,
    uniprot_acc TEXT
);
"""

def extract_uniprot(ted_id: str) -> str:
    # Matches 6 to 10 character UniProt accession (uppercase letters/digits, excluding 'O' and 'Q')
    match = re.search(r"[A-NR-Z0-9]{6,10}", ted_id)
    return match.group(0) if match else ""

def get_db_path(tsv_path: str, db_path: Optional[str] = None) -> Path:
    tsv_path = Path(tsv_path)
    return Path(db_path) if db_path else tsv_path.with_suffix(".duckdb")

def db_exists(db_path: Path) -> bool:
    return db_path.exists()

def create_db(tsv_gz_path: str, db_path: str = None, overwrite: bool = False) -> Path:
    tsv_gz_path = str(tsv_gz_path)
    db_path = str(db_path or Path(tsv_gz_path).with_suffix(".duckdb"))

    if os.path.exists(db_path) and not overwrite:
        print(f"Using existing DB at {db_path}")
        return Path(db_path)

    if os.path.exists(db_path):
        os.remove(db_path)

    con = duckdb.connect(db_path)

    columns = [
        "ted_id", "md5_domain", "consensus_level", "chopping", "nres_domain", "num_segments",
        "plddt", "num_helix_strand_turn", "num_helix", "num_strand", "num_helix_strand", "num_turn",
        "proteome-id", "cath_label", "cath_assignment_level", "cath_assignment_method",
        "packing_density", "norm_rg", "tax_common_name", "tax_scientific_name", "tax_lineage"
    ]

    column_types = {col: "TEXT" for col in columns}

    # Step 1: Read file into DuckDB (streaming, no full DataFrame in memory)
    con.execute("""
        CREATE OR REPLACE TABLE domain_summary AS
        SELECT 
            *,
            regexp_extract("ted_id", 'AF-([A-Z0-9]+)', 1) AS uniprot_acc
        FROM read_csv(
            ?, 
            delim='\t', 
            header=True, 
            columns=?,
            compression='gzip'
        )
    """, [tsv_gz_path, column_types])
    
    con.close()
    print(f"DB created at {db_path}")
    return Path(db_path)

def create_db_slim_copy(tsv_gz_path: str, db_path: str = None, overwrite: bool = False) -> Path:
    tsv_gz_path = str(tsv_gz_path)
    db_path = str(db_path or Path(tsv_gz_path).with_suffix(".slim.duckdb"))

    # ... (same overwrite logic as before) ...
    if os.path.exists(db_path): os.remove(db_path)

    con = duckdb.connect(db_path)
    con.execute("PRAGMA memory_limit='4GB';") # Optional: safety net
    
    # Step 1: Create a table with ONLY the columns you want to keep
    con.execute("""
        CREATE OR REPLACE TABLE domain_summary(
            ted_id VARCHAR,
            chopping VARCHAR,
            cath_label VARCHAR,
            cath_assignment_level VARCHAR,
            cath_assignment_method VARCHAR,
            tax_common_name VARCHAR,
            tax_scientific_name VARCHAR
        );
    """)

    # Step 2: Use COPY to load ONLY those columns from the source file
    # This is highly efficient as it skips reading and processing unused columns.
    path = Path(tsv_gz_path).as_posix()
    con.execute(f"""
        COPY domain_summary FROM '{path}' (
            FORMAT CSV,
            HEADER,
            DELIMITER '\t',
            COMPRESSION GZIP
        );
    """)

    # Step 3: Add and populate the uniprot_acc column
    con.execute("ALTER TABLE domain_summary ADD COLUMN uniprot_acc VARCHAR;")
    con.execute(f"""
        UPDATE domain_summary
        SET uniprot_acc = regexp_extract(ted_id, 'AF-([A-Z0-9]+)', 1);
    """)
    
    con.close()
    print(f"Slim DB created via COPY at {db_path}")
    return Path(db_path)
    
def create_db_from_headerless_tsv(tsv_gz_path: str, db_path: str = None, overwrite: bool = False) -> Path:
    tsv_gz_path = str(tsv_gz_path)
    db_path = str(db_path or Path(tsv_gz_path).with_suffix(".duckdb"))

    if os.path.exists(db_path) and not overwrite:
        print(f"Using existing DB at {db_path}")
        return Path(db_path)
        
    if os.path.exists(db_path):
        os.remove(db_path)

    con = duckdb.connect(db_path)

    # Step 1: Define ALL column names for the source file IN ORDER
    # This is crucial for correctly parsing the header-less file.
    all_column_names = [
        "ted_id", "md5_domain", "consensus_level", "chopping", "nres_domain", "num_segments",
        "plddt", "num_helix_strand_turn", "num_helix", "num_strand", "num_helix_strand", "num_turn",
        "proteome-id", "cath_label", "cath_assignment_level", "cath_assignment_method",
        "packing_density", "norm_rg", "tax_common_name", "tax_scientific_name", "tax_lineage"
    ]
    
    # Define the data types for the columns you are loading.
    # This is still a best practice for performance.
    full_column_types = {
        "ted_id": "VARCHAR", "md5_domain": "VARCHAR", "consensus_level": "VARCHAR",
        "chopping": "VARCHAR", "nres_domain": "INTEGER", "num_segments": "INTEGER",
        "plddt": "FLOAT", "num_helix_strand_turn": "INTEGER", "num_helix": "INTEGER",
        "num_strand": "INTEGER", "num_helix_strand": "INTEGER", "num_turn": "INTEGER",
        "proteome-id": "VARCHAR", "cath_label": "VARCHAR", "cath_assignment_level": "VARCHAR",
        "cath_assignment_method": "VARCHAR", "packing_density": "FLOAT", "norm_rg": "FLOAT",
        "tax_common_name": "VARCHAR", "tax_scientific_name": "VARCHAR", "tax_lineage": "VARCHAR"
    }

    # Step 2: Use read_csv with header=False and provide the names list.
    # Then select the subset of columns you want for your final table.
    con.execute("""
        CREATE OR REPLACE TABLE domain_summary AS
        SELECT 
            -- Explicitly list the columns to keep
            ted_id,
            chopping,
            cath_label,
            cath_assignment_level,
            cath_assignment_method,
            tax_common_name,
            tax_scientific_name,
            -- Also create the new uniprot_acc column
            regexp_extract(ted_id, 'AF-([A-Z0-9]+)', 1) AS uniprot_acc
        FROM read_csv(
            ?, 
            delim='\t', 
            header=False,         -- CRITICAL: Set header to False
            columns=?,            -- Provide types for efficiency
            compression='gzip'
        )
    """, [tsv_gz_path, full_column_types])
    
    con.close()
    print(f"Slim DB created from header-less TSV at {db_path}")
    return Path(db_path)


def query_by_uniprot_ids(db_path: str, uniprot_ids: List[str], keyword: Optional[str] = None) -> pd.DataFrame:
    con = duckdb.connect(db_path)
    id_list_str = ', '.join(f"'{id}'" for id in uniprot_ids)
    query = f"""
        SELECT * FROM domain_summary
        WHERE uniprot_acc IN ({id_list_str})
    """
    if keyword:
        query += f" AND tax_common_name ILIKE '%{keyword}%'"

    df = con.execute(query).fetchdf()
    con.close()
    return df