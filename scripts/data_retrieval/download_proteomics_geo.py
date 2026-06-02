import os
import urllib.request
import urllib.error
import ssl
import sys
import tarfile
import json
import re

def download_gse_matrix(gse_id, outdir):
    print(f"Downloading supplementary matrix for {gse_id}...")
    os.makedirs(outdir, exist_ok=True)
    
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    
    # Use Entrez to find the exact FTP URL
    try:
        gse_num = gse_id.replace("GSE", "")
        uid_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse_id}[Accession]&retmode=json"
        
        req = urllib.request.Request(uid_url)
        with urllib.request.urlopen(req, context=ctx) as r:
            uid_data = json.loads(r.read())
            uid = uid_data['esearchresult']['idlist'][0]
            
        sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={uid}&retmode=json"
        req = urllib.request.Request(sum_url)
        
        with urllib.request.urlopen(req, context=ctx) as r:
            sum_data = json.loads(r.read())
            ftp_base = sum_data['result'][uid]['ftplink']
            
            # GEO supplementary files are always under supp/
            url = f"{ftp_base}suppl/{gse_id}_RAW.tar"
            print(f"Resolved FTP Path: {url}")
            
    except Exception as e:
        print(f"Failed to resolve FTP path: {e}")
        return
    
    try:
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0')
        req.add_header('Accept', '*/*')
        
        dest_tar = os.path.join(outdir, f"{gse_id}_RAW.tar")
        
        print(f"Fetching from {url}...")
        with urllib.request.urlopen(req, context=ctx) as response, open(dest_tar, 'wb') as out_file:
            # chunked download
            chunk_size = 1024 * 1024 * 10 # 10MB
            total_size = 0
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
                total_size += len(chunk)
                sys.stdout.write(f"\rDownloaded: {total_size / (1024*1024):.1f} MB")
                sys.stdout.flush()
        print(f"\nSuccessfully downloaded to {dest_tar}")
        
    except Exception as e:
        print(f"\nFailed to download {gse_id}: {e}")

if __name__ == "__main__":
    outdir = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/Proteomics/GSE276114"
    download_gse_matrix("GSE276114", outdir)
    print("\nExtraction of RAW tarball is required post-download.")
