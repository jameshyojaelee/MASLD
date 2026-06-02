import os
import sys
import argparse
import requests

def fetch_metadata(gse_id, out_dir):
    print(f"Fetching metadata for {gse_id}...")
    
    # Files
    run_info_file = os.path.join(out_dir, f"{gse_id}_run_info.csv")
    
    # SRA Run Info URL
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={gse_id}"
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    
    try:
        print(f"Downloading RunInfo from: {url}")
        response = requests.get(url, headers=headers, allow_redirects=True)
        response.raise_for_status()
        
        content = response.text
        
        # Check if valid CSV header
        if "Run,ReleaseDate" not in content and "Run" not in content:
            print(f"Error: Unexpected response format:\n{content[:200]}")
            sys.exit(1)
            
        with open(run_info_file, "w") as f:
            f.write(content)
            
        print(f"Saved {run_info_file}")
        
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gse_id")
    parser.add_argument("out_dir")
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    fetch_metadata(args.gse_id, args.out_dir)
