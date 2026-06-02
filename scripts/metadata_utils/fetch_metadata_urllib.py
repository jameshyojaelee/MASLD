import urllib.request
import urllib.parse
import urllib.error
import json
import time
import os
import csv
import sys
import re

# Constants
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GSE_LIST = ['GSE162876', 'GSE213621', 'GSE274914', 'GSE224069', 'GSE263273', 'GSE225616']
OUTPUT_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"

def get_json(url):
    time.sleep(0.5) # Rate limit
    try:
        with urllib.request.urlopen(url) as response:
            return json.loads(response.read().decode())
    except Exception as e:
        print(f"Error fetching JSON from {url}: {e}")
        return None

def esearch(db, term):
    params = {
        'db': db,
        'term': term,
        'retmode': 'json',
        'retmax': 1000
    }
    url = BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode(params)
    data = get_json(url)
    if data:
        return data.get('esearchresult', {}).get('idlist', [])
    return []

def esummary(db, id_list):
    ids_str = ",".join(id_list)
    params = {
        'db': db,
        'id': ids_str,
        'retmode': 'json'
    }
    url = BASE_URL + "esummary.fcgi?" + urllib.parse.urlencode(params)
    return get_json(url)

def efetch_runinfo(ids):
    if not ids:
        return None
    ids_str = ",".join(ids)
    params = {
        'db': 'sra',
        'id': ids_str,
        'retmode': 'text',
        'rettype': 'runinfo'
    }
    url = BASE_URL + "efetch.fcgi?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url) as response:
            return response.read().decode()
    except Exception as e:
        print(f"Error fetching runinfo: {e}")
        return None

def main():
    print(f"Starting metadata fetch for {len(GSE_LIST)} datasets...")
    
    total_samples = 0
    
    for gse in GSE_LIST:
        print(f"\nProcessing {gse}...")
        gse_dir = os.path.join(OUTPUT_DIR, gse)
        meta_dir = os.path.join(gse_dir, "metadata")
        os.makedirs(meta_dir, exist_ok=True)
        
        # 1. Search GDS for GSE
        gds_ids = esearch('gds', f"{gse}[Accession]")
        if not gds_ids:
            print(f"  GSE not found in GDS: {gse}")
            continue
            
        # 2. Get Summary to find SRP
        summary_json = esummary('gds', gds_ids)
        srp_id = None
        if summary_json and 'result' in summary_json:
            for uid in gds_ids:
                if uid in summary_json['result']:
                    item = summary_json['result'][uid]
                    # Look for ExtRelation or typical text fields
                    # The summary structure varies, but often contains 'extrelations' list
                    if 'extrelations' in item:
                        for rel in item['extrelations']:
                            if rel.get('relationtype') == 'SRA':
                                srp_id = rel.get('targetobject')
                                break
        
        if not srp_id:
             print(f"  Could not find SRP/SRA link for {gse} in GDS summary.")
             # Fallback: try searching SRA with the GSE accession as term again (sometimes works differently)
             # But we verified it failed.
             continue
             
        print(f"  Found SRA Study: {srp_id}")
        
        # 3. Search SRA for SRP
        sra_ids = esearch('sra', srp_id)
        print(f"  Found {len(sra_ids)} SRA entries for {srp_id}")
        
        # 4. Fetch RunInfo
        runinfo_csv = efetch_runinfo(sra_ids)
        
        if runinfo_csv:
            filename = os.path.join(meta_dir, f"{gse}_run_info.csv")
            with open(filename, 'w') as f:
                f.write(runinfo_csv)
            
            # Count runs (skip header)
            lines = runinfo_csv.strip().split('\n')
            run_count = len(lines) - 1 if len(lines) > 0 else 0
            print(f"  Saved run info to {filename} ({run_count} runs)")
            total_samples += run_count
        else:
            print("  Failed to retrieve run info.")
            
    print(f"\nTotal samples found: {total_samples}")

if __name__ == "__main__":
    main()
