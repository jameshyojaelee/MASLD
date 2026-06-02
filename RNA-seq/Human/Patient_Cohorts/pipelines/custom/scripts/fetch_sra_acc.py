import requests
import xml.etree.ElementTree as ET
import sys
import time

def get_sra_ids(term):
    # Search for the term
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "sra",
        "term": term,
        "retmax": 500,
        "usehistory": "y"
    }
    response = requests.get(search_url, params=params)
    response.raise_for_status()
    
    root = ET.fromstring(response.content)
    count = int(root.find("Count").text)
    webenv = root.find("WebEnv").text
    query_key = root.find("QueryKey").text
    
    print(f"Found {count} entries for {term}")
    
    # Fetch run info
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {
        "db": "sra",
        "WebEnv": webenv,
        "query_key": query_key,
        "retmode": "xml"
    }
    
    # Fetching in batches might be safer but for <500 one shot is fine
    response = requests.get(fetch_url, params=fetch_params)
    response.raise_for_status()
    
    # Parse the XML to extract run accessions
    # The output structure is complex, containing Experiment, Run, etc.
    # We want info from <RUN accession="...">
    # Check if we can get runinfo text format directly
    
    fetch_params['rettype'] = 'runinfo'
    fetch_params['retmode'] = 'text'
    
    response = requests.get(fetch_url, params=fetch_params)
    response.raise_for_status()
    
    lines = response.text.strip().split('\n')
    header = lines[0].split(',')
    
    try:
        run_idx = header.index("Run")
    except ValueError:
        print("Could not find 'Run' in runinfo header")
        return []

    run_accessions = []
    for line in lines[1:]:
        parts = line.split(',')
        if len(parts) > run_idx:
            run_accessions.append(parts[run_idx])
            
    return run_accessions

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fetch_sra_acc.py <SEARCH_TERM> [OUTPUT_FILE]")
        sys.exit(1)
        
    term = sys.argv[1]
    
    try:
        accs = get_sra_ids(term)
        if len(sys.argv) >= 3:
            with open(sys.argv[2], 'w') as f:
                for acc in accs:
                    f.write(f"{acc}\n")
            print(f"Wrote {len(accs)} accessions to {sys.argv[2]}")
        else:
            for acc in accs:
                print(acc)
    except Exception as e:
        print(f"Error fetching data: {e}")
        sys.exit(1)
