import requests
import time

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