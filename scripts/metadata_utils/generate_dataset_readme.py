import os
import csv
import glob
import urllib.request
import urllib.parse
import json
import time

BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
OUTPUT_README = os.path.join(BASE_DIR, "README.md")
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get_json(url):
    time.sleep(0.34)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            return json.loads(response.read().decode())
    except Exception as e:
        print(f"Error fetching JSON from {url}: {e}")
        return None

def fetch_gse_summary(gse_id):
    # Search GDS to get Summary
    search_url = BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode({'db': 'gds', 'term': f"{gse_id}[Accession]", 'retmode': 'json'})
    search_data = get_json(search_url)
    if not search_data or 'esearchresult' not in search_data or not search_data['esearchresult']['idlist']:
        return "No summary available", "No title available"
    
    uid = search_data['esearchresult']['idlist'][0]
    summary_url = BASE_URL + "esummary.fcgi?" + urllib.parse.urlencode({'db': 'gds', 'id': uid, 'retmode': 'json'})
    summary_data = get_json(summary_url)
    
    if summary_data and 'result' in summary_data and uid in summary_data['result']:
        res = summary_data['result'][uid]
        title = res.get('title', 'No title')
        summary = res.get('summary', 'No summary')
        return title, summary
    return "No summary available", "No title available"

def get_stats(gse_path, gse_id):
    run_info = glob.glob(os.path.join(gse_path, "metadata", "*_run_info.csv"))
    attr_info = glob.glob(os.path.join(gse_path, "metadata", "*_sample_attributes.csv"))
    
    n_runs = 0
    details = {
        'Diets': set(),
        'Tissues': set(),
        'Strains': set(),
        'Genotypes': set(),
        'DirectTreatments': set(),
        'Timepoints': set()
    }
    
    if run_info:
        with open(run_info[0]) as f:
            n_runs = len(f.readlines()) - 1
            
    if attr_info:
        with open(attr_info[0]) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Normalize values
                for k, v in row.items():
                    val_lower = v.lower()
                    if k in ['diet', 'treatment'] or 'diet' in val_lower or 'hfd' in val_lower or 'chow' in val_lower:
                        if 'hfd' in val_lower or 'high fat' in val_lower: details['Diets'].add(v)
                        elif 'chow' in val_lower: details['Diets'].add(v)
                        elif 'amln' in val_lower: details['Diets'].add(v)
                        elif 'western' in val_lower: details['Diets'].add(v)
                        elif 'control' in val_lower: details['Diets'].add(v)
                        else: details['DirectTreatments'].add(v)
                        
                    if k in ['tissue', 'source_name']:
                        details['Tissues'].add(v)
                    if k in ['strain', 'organism']:
                        details['Strains'].add(v)
                    if k in ['genotype']:
                        details['Genotypes'].add(v)
                    if k in ['timepoint', 'time point', 'age', 'duration']:
                        details['Timepoints'].add(v)

    # Fetch Title/Summary
    title, summary = fetch_gse_summary(gse_id)
    
    return n_runs, details, title, summary

def main():
    datasets = []
    
    print("Gathering info...")
    for gse in sorted(os.listdir(BASE_DIR)):
        gse_path = os.path.join(BASE_DIR, gse)
        if not os.path.isdir(gse_path) or not gse.startswith("GSE"): continue
        
        print(f"  Processing {gse}...")
        n, details, title, summary = get_stats(gse_path, gse)
        datasets.append({
            'id': gse,
            'runs': n,
            'details': details,
            'title': title,
            'summary': summary
        })

    print("Generating README...")
    with open(OUTPUT_README, 'w') as f:
        f.write("# Other Diet Mouse RNA-seq Datasets\n\n")
        f.write("Repository of external mouse RNA-seq datasets for MASLD/NASH research, specifically focusing on Diet models (HFD, Western, AMLN) excluding Choline Deficiency (MCD).\n\n")
        f.write("## Overview\n\n")
        f.write("| Dataset | Runs | Diet/Model | Tissue | Summary |\n")
        f.write("|---|---|---|---|---|\n")
        
        for d in datasets:
            diets = ", ".join(sorted(d['details']['Diets'])) if d['details']['Diets'] else "Unspecified"
            tissues = ", ".join(sorted(d['details']['Tissues'])) if d['details']['Tissues'] else "Unspecified"
            # Truncate summary for table
            short_sum = (d['title'][:100] + '...') if len(d['title']) > 100 else d['title']
            f.write(f"| [{d['id']}](./{d['id']}) | {d['runs']} | {diets} | {tissues} | {short_sum} |\n")
            
        f.write("\n## Detailed Dataset Metadata\n")
        for d in datasets:
            f.write(f"\n### [{d['id']}](./{d['id']})\n")
            f.write(f"**Title**: {d['title']}\n\n")
            f.write(f"**Runs**: {d['runs']}\n\n")
            
            f.write("**Key Attributes**:\n")
            det = d['details']
            if det['Diets']: f.write(f"- **Diets**: {', '.join(sorted(det['Diets']))}\n")
            if det['Tissues']: f.write(f"- **Tissues**: {', '.join(sorted(det['Tissues']))}\n")
            if det['Strains']: f.write(f"- **Strains**: {', '.join(sorted(det['Strains']))}\n")
            if det['Genotypes']: f.write(f"- **Genotypes**: {', '.join(sorted(det['Genotypes']))}\n")
            if det['Timepoints']: f.write(f"- **Timepoints/Age**: {', '.join(sorted(det['Timepoints']))}\n")
            if det['DirectTreatments']: f.write(f"- **Treatments**: {', '.join(sorted(det['DirectTreatments']))}\n")
            
            f.write(f"\n**Summary**:\n{d['summary']}\n")
            f.write("\n---\n")

    print(f"Generated {OUTPUT_README}")

if __name__ == "__main__":
    main()
