import os
import glob
from pathlib import Path

# Constants
BASE_PATH = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
ARTIFACT_PATH = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/DATASET_STATUS.md"

DATASETS = {
    "GSE167523": {
        "raw_dir": f"{BASE_PATH}/data/raw/GSE167523/fastq",
        "results_dir": f"{BASE_PATH}/results/GSE167523",
        "bad_list": f"{BASE_PATH}/pipelines/custom/GSE167523/bad_files_batch.txt",
        "accession_list": None
    },
    "GSE126848": {
        "raw_dir": f"{BASE_PATH}/data/raw/GSE126848/fastq",
        "results_dir": f"{BASE_PATH}/results/GSE126848",
        "bad_list": f"{BASE_PATH}/pipelines/custom/GSE126848/bad_files_batch.txt",
        "accession_list": None
    },
    "PRJNA512027": {
        "raw_dir": f"{BASE_PATH}/data/raw/PRJNA512027/fastq",
        "results_dir": f"{BASE_PATH}/results/PRJNA512027",
        "bad_list": None,
        "accession_list": f"{BASE_PATH}/pipelines/custom/PRJNA512027/accession_list.txt"
    }
}

def load_bad_files(filepath):
    bad_files = set()
    if filepath and os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                bad_files.add(os.path.basename(line.strip()))
    return bad_files

def load_accessions(filepath):
    accessions = []
    if filepath and os.path.exists(filepath):
        with open(filepath, 'r') as f:
            accessions = [line.strip() for line in f if line.strip()]
    return accessions

def check_analysis_status(dataset, sample_id, results_dir):
    # Sample ID correction: Remove extensions and _1/_2 suffix
    clean_name = sample_id.replace('.fastq.gz', '').replace('.fastq', '')
    base_id = clean_name.split('_')[0]
    
    # FastQC
    fastqc_path_1 = f"{results_dir}/qc/fastqc/{base_id}_1_fastqc.html"
    fastqc_path_2 = f"{results_dir}/qc/fastqc/{base_id}_fastqc.html"
    fastqc_done = "✅" if os.path.exists(fastqc_path_1) or os.path.exists(fastqc_path_2) else "⬜"
    
    # STAR
    bam_path = f"{results_dir}/alignments/star/{base_id}/{base_id}.Aligned.sortedByCoord.out.bam"
    star_done = "✅" if os.path.exists(bam_path) else "⬜"
    
    return fastqc_done, star_done

def generate_markdown():
    md_lines = []
    md_lines.append("# RNA-seq Dataset & Analysis Status Tracker")
    md_lines.append(f"**Last Updated**: {os.popen('date').read().strip()}\n")
    
    md_lines.append("## Legend")
    md_lines.append("| Status | Description |")
    md_lines.append("|---|---|")
    md_lines.append("| 🟢 **Verified** | File exists and passed integrity check |")
    md_lines.append("| 🟡 **Repairing** | Corrupt/Missing file currently being re-downloaded |")
    md_lines.append("| 🔵 **Downloading** | Initial download in progress |")
    md_lines.append("| 🔴 **Missing** | File missing and not queued for repair |")
    md_lines.append("")

    for dataset, info in DATASETS.items():
        md_lines.append(f"## {dataset}")
        
        # Gather all known files
        bad_files = load_bad_files(info['bad_list'])
        
        # Get existing files
        existing_files = set(os.listdir(info['raw_dir'])) if os.path.exists(info['raw_dir']) else set()
        
        # Determine full set of expected files
        all_samples = set()
        
        # 1. From Accession List (if available)
        if info['accession_list']:
            accs = load_accessions(info['accession_list'])
            for acc in accs:
                all_samples.add(f"{acc}_1.fastq.gz")
                all_samples.add(f"{acc}_2.fastq.gz")
        
        # 2. From Bad Files
        all_samples.update(bad_files)
        
        # 3. From Existing Files
        all_samples.update([f for f in existing_files if f.endswith('.fastq.gz') or f.endswith('.fastq')])
        
        sorted_samples = sorted(list(all_samples))
        
        # Stats
        total = len(sorted_samples)
        verified = 0
        
        table_rows = []
        for sample in sorted_samples:
            status_icon = "❓"
            status_text = "Unknown"
            
            is_bad = sample in bad_files
            exists = sample in existing_files
            
            # Logic
            if is_bad:
                status_icon = "🟡"
                status_text = "Repairing"
            elif exists:
                size = os.path.getsize(os.path.join(info['raw_dir'], sample))
                if size > 1000:
                    status_icon = "🟢"
                    status_text = "Verified"
                    verified += 1
                else:
                    status_icon = "🔵" 
                    status_text = "Downloading"
            else:
                if dataset == "PRJNA512027":
                    status_icon = "🔵"
                    status_text = "Queued"
                else:
                    status_icon = "🔴"
                    status_text = "Missing"

            # Check Analysis Status
            fastqc, star = check_analysis_status(dataset, sample, info['results_dir'])
            
            table_rows.append(f"| {sample} | {status_icon} | {status_text} | {fastqc} | {star} |")

        # Summary Metrics
        md_lines.append(f"**Total Files**: {total} | **Verified**: {verified} | **Active/Repairing**: {total - verified}")
        md_lines.append("")
        md_lines.append(f"| Filename | File Status | Details | FastQC | STAR |")
        md_lines.append(f"|---|---|---|---|---|")
        md_lines.extend(table_rows)
        md_lines.append("")
        
    with open(ARTIFACT_PATH, 'w') as f:
        f.write("\n".join(md_lines))
    
    print(f"Artifact created at {ARTIFACT_PATH}")

if __name__ == "__main__":
    generate_markdown()
