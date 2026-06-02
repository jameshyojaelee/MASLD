import os
import glob
import csv

BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
JOB_LIST = os.path.join(BASE_DIR, "master_download_list.txt")

def main():
    jobs = []
    
    # scan for run_info.csv
    # e.g. .../GSE12345/metadata/GSE12345_run_info.csv
    for gse_dir in glob.glob(os.path.join(BASE_DIR, "GSE*")):
        gse_id = os.path.basename(gse_dir)
        meta = os.path.join(gse_dir, "metadata", f"{gse_id}_run_info.csv")
        
        if os.path.exists(meta):
            with open(meta, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'Run' in row:
                        jobs.append(f"{gse_id}\t{row['Run']}")
    
    with open(JOB_LIST, 'w') as f:
        f.write("\n".join(jobs))
        f.write("\n")
        
    print(f"Generated {len(jobs)} jobs in {JOB_LIST}")

if __name__ == "__main__":
    main()
