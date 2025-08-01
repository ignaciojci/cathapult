import requests
import time
import gzip
import pandas as pd
from io import StringIO
from tqdm import tqdm


def fetch_ted_summary(uniprot_id, timeout=10, delay=1):
    url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}?skip=0&limit=100"
    try:
        response = requests.get(url, headers={"accept": "application/json"}, timeout=timeout)
        response.raise_for_status()
        data = response.json().get("data", [])
        time.sleep(delay)  # optional delay
        return data
    except requests.exceptions.Timeout:
        print(f"Timeout for {uniprot_id}")
    except requests.exceptions.RequestException as e:
        print(f"Request failed for {uniprot_id}: {e}")
    return []

def get_ted_domains():
    path = os.getenv("TED_DOMAINS")
    if path is None:
        raise EnvironmentError(
            "Environment variable 'TED_DOMAINS' not set.\n"
            "Please set it to the path of 'ted_365m.domain_summary.cath.globularity.taxid.tsv.gz'."
        )
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    return path
    

def extract_uniprot_id(model_id):
    """Extract UniProt ID from something like 'AF-A0A000-F1-model_v4_TED01'."""
    parts = model_id.split('-')
    if len(parts) >= 2:
        return parts[1]
    return None
    
def open_maybe_gzip(path):
    with open(path, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':  # gzip magic number
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')

TED_COLUMNS = [
    "ted_id", "md5_domain", "consensus_level", "chopping", "nres_domain",
    "num_segments", "plddt", "num_helix_strand_turn", "num_helix", "num_strand",
    "num_helix_strand", "num_turn", "proteome-id", "cath_label",
    "cath_assignment_level", "cath_assignment_method", "packing_density",
    "norm_rg", "tax_common_name", "tax_scientific_name", "tax_lineage"
]

def filter_ted_summary(
    target_ids,
    ted_db_gz: str,
    filter_keyword: str = "sapiens",
    total_lines: int = 364_806_077
) -> pd.DataFrame:
    # Collect matching lines in memory
    lines = []
    output_ids = []
    with open_maybe_gzip(ted_db_gz) as fin:
        for line in tqdm(fin, total=total_lines, unit='lines', ncols=100, desc="Filtering", smoothing=0.01):
            if filter_keyword not in line:
                continue
            col1 = line.split('\t', 1)[0]
            uniprot_id = extract_uniprot_id(col1)
            if uniprot_id in target_ids:
                lines.append(line)
                output_ids.append(uniprot_id)

    # Combine into one string and load into pandas
    if lines:
        data_str = ''.join(lines)
        df = pd.read_csv(StringIO(data_str), sep='\t', names=TED_COLUMNS)
        df['uniprot_id'] = output_ids
        return df
    else:
        return pd.DataFrame(columns=TED_COLUMNS)  # Empty if no match