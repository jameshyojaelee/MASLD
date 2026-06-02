import os
import ssl
import sys
import time
import argparse
import urllib.request
import urllib.parse
import json

def fetch_with_retry(url, retries=5):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    for i in range(retries):
        try:
            time.sleep(1) # Honor rate limit
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, context=ctx, timeout=30) as res:
                return res.read().decode('utf-8')
        except Exception as e:
            print(f"Warning: Attempt {i+1} failed: {e}")
            time.sleep(2 ** i)
    raise Exception(f"Failed to fetch {url} after {retries} attempts")

def get_srr_for_gse(gse, outdir):
    print(f"Fetching SRR IDs for {gse}...")
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # 0. Convert GSE to GDS UID
    try:
        gse_num = int(gse.replace("GSE", ""))
        gds_uid = 200000000 + gse_num
    except ValueError:
        print(f"Invalid GSE format: {gse}")
        return

    # 1. Fetch SRA accession (SRP...) from GDS
    summary_url = f"{base_url}esummary.fcgi?db=gds&id={gds_uid}&retmode=json"
    summary_data = json.loads(fetch_with_retry(summary_url))
    
    srp_acc = None
    try:
        extrels = summary_data["result"][str(gds_uid)].get("extrelations", [])
        for rel in extrels:
            if rel.get("relationtype") == "SRA":
                srp_acc = rel.get("targetobject")
                break
    except KeyError:
        pass

    if not srp_acc:
        print(f"Warning: No SRA relation found for {gse}.")
        return
        
    print(f"Mapped {gse} -> {srp_acc}. Querying SRA...")

    # 2. esearch in sra database using SRP
    search_url = f"{base_url}esearch.fcgi?db=sra&term={srp_acc}[Accession]&retmax=5000&retmode=json"
    search_data = json.loads(fetch_with_retry(search_url))
    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    
    if not id_list:
        print(f"No SRA records found for {srp_acc}.")
        return

    # 3. efetch runinfo (in batches to avoid URI too long)
    batch_size = 200
    srr_list = set()
    
    for i in range(0, len(id_list), batch_size):
        batch_ids = ",".join(id_list[i:i+batch_size])
        fetch_url = f"{base_url}efetch.fcgi?db=sra&id={batch_ids}&rettype=runinfo&retmode=csv"
        
        csv_data = fetch_with_retry(fetch_url)
        lines = csv_data.strip().split('\n')
        
        for line in lines[1:]: # skip header
            cols = line.split(',')
            if cols and (cols[0].startswith('SRR') or cols[0].startswith('ERR')):
                srr_list.add(cols[0])
                
    if srr_list:
        os.makedirs(outdir, exist_ok=True)
        out_file = os.path.join(outdir, f"{gse}_SRR_list.txt")
        with open(out_file, 'w') as f:
            for srr in sorted(list(srr_list)):
                f.write(f"{srr}\n")
        print(f"Saved {len(srr_list)} SRRs to {out_file}")
    else:
        print(f"Warning: Found SRA records but no SRR IDs for {gse}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('gse_list', nargs='+')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    
    for gse in args.gse_list:
        try:
            get_srr_for_gse(gse, args.outdir)
        except Exception as e:
            print(f"Error processing {gse}: {e}")
