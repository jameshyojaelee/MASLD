import urllib.request
import urllib.parse
import urllib.error
import json
import time
import os
import csv
import sys
import xml.etree.ElementTree as ET

# Constants
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GSE_LIST = ['GSE162876', 'GSE213621', 'GSE274914', 'GSE224069', 'GSE263273', 'GSE225616']
OUTPUT_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"

def get_json(url):
    for i in range(3):
        time.sleep(1 + i) # Increasing backoff
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                return json.loads(response.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                print(f"  Rate limited (429). Retrying {i+1}/3...")
                time.sleep(2 + i*2)
                continue
            print(f"Error fetching JSON from {url}: {e}")
            return None
        except Exception as e:
            print(f"Error fetching JSON from {url}: {e}")
            return None
    return None

# ... (omitted xml/search/summary functions if unchanged) ...

def main():
    print(f"Starting detailed metadata fetch for {len(GSE_LIST)} datasets...")
    
    total_runs = 0
    
    for gse in GSE_LIST:
        print(f"Processing {gse}...")
        gse_dir = os.path.join(OUTPUT_DIR, gse)
        meta_dir = os.path.join(gse_dir, "metadata")
        os.makedirs(meta_dir, exist_ok=True)
        
        # Strategy: GSE -> GDS (Summary) -> SRP -> SRA (RunInfo) -> BioProject/BioSample (Attributes)
        
        # 1. GSE -> SRP
        gds_ids = esearch('gds', f"{gse}[Accession]")
        srp_id = None
        if gds_ids:
            # OPTIMIZATION: Only check first 5 GDS entries to find the linked Study/SRP.
            # Checking all hundreds leads to 429 errors.
            subset_ids = gds_ids[:5]
            summary = esummary('gds', subset_ids)

def get_xml(url):
    time.sleep(0.34)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            return ET.fromstring(response.read().decode())
    except Exception as e:
        print(f"Error fetching XML from {url}: {e}")
        return None

def esearch(db, term, retmax=1000):
    params = {
        'db': db,
        'term': term,
        'retmode': 'json',
        'retmax': retmax
    }
    url = BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode(params)
    data = get_json(url)
    if data:
        return data.get('esearchresult', {}).get('idlist', [])
    return []

def esummary(db, id_list):
    if not id_list:
        return {}
    ids_str = ",".join(id_list)
    params = {
        'db': db,
        'id': ids_str,
        'retmode': 'json'
    }
    url = BASE_URL + "esummary.fcgi?" + urllib.parse.urlencode(params)
    data = get_json(url)
    if data and 'result' in data:
        return data['result']
    return {}

def efetch(db, ids):
    if not ids:
        return None
    ids_str = ",".join(ids)
    # Using XML for detailed SRA metadata (sample attributes usually in XML)
    params = {
        'db': db,
        'id': ids_str,
        'retmode': 'xml' 
    }
    url = BASE_URL + "efetch.fcgi?" + urllib.parse.urlencode(params)
    return get_xml(url)

def fetch_sra_metadata_flat(sra_ids):
    """
    Fetches runinfo (flat CSV) for a list of SRA/SRP IDs. (Efficient for SRR info)
    """
    ids_str = ",".join(sra_ids)
    params = {
        'db': 'sra',
        'id': ids_str,
        'retmode': 'text',
        'rettype': 'runinfo'
    }
    url = BASE_URL + "efetch.fcgi?" + urllib.parse.urlencode(params)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            return response.read().decode()
    except Exception as e:
        print(f"Error fetching runinfo: {e}")
        return None

def parse_biosample_attrs(biosample_id):
    """
    Fetches BioSample text to parse attributes like 'diet', 'tissue', 'treatment'
    """
    xml = efetch('biosample', [biosample_id])
    attrs = {}
    if xml:
        for attr in xml.findall('.//Attribute'):
            key = attr.attrib.get('attribute_name', '').lower()
            val = attr.text
            if key and val:
                attrs[key] = val
    return attrs

def main():
    print(f"Starting detailed metadata fetch for {len(GSE_LIST)} datasets...")
    
    total_runs = 0
    
    for gse in GSE_LIST:
        print(f"Processing {gse}...")
        gse_dir = os.path.join(OUTPUT_DIR, gse)
        meta_dir = os.path.join(gse_dir, "metadata")
        os.makedirs(meta_dir, exist_ok=True)
        
        # Strategy: GSE -> GDS (Summary) -> SRP -> SRA (RunInfo) -> BioProject/BioSample (Attributes)
        
        # 1. GSE -> SRP
        gds_ids = esearch('gds', f"{gse}[Accession]")
        srp_id = None
        if gds_ids:
            summary = esummary('gds', gds_ids)
            for uid in gds_ids:
                if uid in summary:
                    item = summary[uid]
                    if 'extrelations' in item:
                        for rel in item['extrelations']:
                            if rel.get('relationtype') == 'SRA':
                                srp_id = rel.get('targetobject') 
                                break # Found one SRP
                    if srp_id: break
        
        # FIX: If srp_id is actually an SRX (Experiment), we need to find the parent SRP (Study)
        if srp_id and srp_id.startswith('SRX'):
            print(f"  {srp_id} is an Experiment, not a Study. Resolving parent SRP...")
            # Fetch SRX summary or runinfo to find SRP
            srx_info = fetch_sra_metadata_flat([srp_id])
            if srx_info:
                reader = csv.DictReader(srx_info.strip().split('\n'))
                rows = list(reader)
                if rows and 'SRAStudy' in rows[0]:
                    srp_id = rows[0]['SRAStudy']
                    print(f"  Resolved parent SRP: {srp_id}")
                    
        if not srp_id:
            # Fallback: Search SRA directly for GSE nickname
            print(f"  GDS Link failed. Searching SRA for {gse} directly...")
            sra_search = esearch('sra', gse)
            if sra_search:
                print(f"  Found {len(sra_search)} SRA entries via direct search.")
                # We need to find the shared Study ID from these entries
                # Fetch RunInfo for the first few to identify the SRAStudy (SRP)
                first_batch = sra_search[:5] if len(sra_search) > 5 else sra_search
                temp_info = fetch_sra_metadata_flat(first_batch)
                if temp_info:
                    reader = csv.DictReader(temp_info.strip().split('\n'))
                    rows = list(reader)
                    if rows and 'SRAStudy' in rows[0]:
                        srp_id = rows[0]['SRAStudy']
                        print(f"  Identified common SRP from direct search: {srp_id}")
            else:
                print(f"  FAILED to find entry for {gse}")
                continue
        else:
             print(f"  Found SRP via GDS: {srp_id}")

        # Now, ensuring we have the SRP, fetch ALL runs for this Study
        if srp_id:
             print(f"  Fetching comprehensive run list for Study {srp_id}...")
             sra_search = esearch('sra', srp_id, retmax=5000)
             run_info = fetch_sra_metadata_flat(sra_search)
        else:
             print(f"  Could not identify SRP for {gse}. Skipping.")
             continue

        if run_info:
            run_csv_path = os.path.join(meta_dir, f"{gse}_run_info.csv")
            with open(run_csv_path, 'w') as f:
                f.write(run_info)
            
            # Now Parse Attributes for README
            # 1. Read CSV
            reader = csv.DictReader(run_info.strip().split('\n'))
            runs = list(reader)
            print(f"  Saved {len(runs)} runs.")
            total_runs += len(runs)
            
            # 2. Get Unique BioSample IDs
            biosamples = list(set([r['BioSample'] for r in runs if 'BioSample' in r]))
            print(f"  Fetching attributes for {len(biosamples)} BioSamples...")
            
            # 3. Fetch BioSample Attributes (Batch)
            # Eutils allows batch fetching.
            # We will generate a 'sample_attributes.csv'
            
            attr_table = []
            
            # Batch in 50
            batch_size = 50
            for i in range(0, len(biosamples), batch_size):
                batch = biosamples[i:i+batch_size]
                # XML fetch
                xml_root = efetch('biosample', batch)
                if xml_root is None: continue
                
                # Parse
                for sample in xml_root.findall('BioSample'):
                    bs_id = sample.attrib.get('accession')
                    row = {'BioSample': bs_id}
                    
                    # Common keys: tissue, diet, treatment, strain, age, sex
                    for attr in sample.findall('.//Attribute'):
                        k = attr.attrib.get('attribute_name', '').lower()
                        v = attr.text
                        if k:
                            row[k] = v
                    attr_table.append(row)
            
            # Save Attribute Table
            if attr_table:
                # Get all headers
                all_keys = set().union(*(d.keys() for d in attr_table))
                # Prioritize important ones
                priority = ['BioSample', 'source_name', 'tissue', 'diet', 'treatment', 'genotype', 'strain', 'age', 'sex']
                fieldnames = [c for c in priority if c in all_keys] + sorted([c for c in all_keys if c not in priority])
                
                attr_csv_path = os.path.join(meta_dir, f"{gse}_sample_attributes.csv")
                with open(attr_csv_path, 'w') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(attr_table)
                print(f"  Saved detailed attributes to {attr_csv_path}")

    print(f"\nTotal Runs Processed: {total_runs}")

if __name__ == "__main__":
    main()
