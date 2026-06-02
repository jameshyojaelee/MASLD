import os
import argparse
import subprocess
import json

def get_srr_metadata(gse_list, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    for gse in gse_list:
        print(f"Fetching metadata for {gse}...")
        try:
            # Use ffq to get metadata for GSE
            result = subprocess.run(['ffq', '--ftp', gse], capture_output=True, text=True, check=True)
            data = json.loads(result.stdout)
            
            srr_list = []
            for item in data:
                # ffq structure can vary, look for SRR runs
                if 'runs' in item:
                    for run_id, run_info in item['runs'].items():
                        srr_list.append(run_id)
                elif 'accession' in item and item['accession'].startswith('SRR'):
                     srr_list.append(item['accession'])
                     
            if not srr_list:
                print(f"Warning: No SRR found for {gse} via ffq.")
            else:
                out_file = os.path.join(output_dir, f"{gse}_SRR_list.txt")
                with open(out_file, 'w') as f:
                    for srr in sorted(list(set(srr_list))):
                        f.write(f"{srr}\n")
                print(f"Saved {len(srr_list)} SRRs to {out_file}")
                
        except Exception as e:
            print(f"Error fetching {gse}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch SRR lists for GSE accessions")
    parser.add_argument('--gse', nargs='+', required=True, help='List of GSE accessions')
    parser.add_argument('--outdir', required=True, help='Output directory for SRR lists')
    args = parser.parse_args()
    
    get_srr_metadata(args.gse, args.outdir)
