import pandas as pd
import os

# Define Paths
BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"

# Map Dataset Name -> Metadata File Path
dataset_paths = {
    # Human
    "Govaere (GSE135251)": "Human/Patient_Cohorts/archive/old_data_metadata/metadata/GSE135251_SraRunTable.csv",
    "Gerhard (PRJNA512027)": "Human/Patient_Cohorts/pipelines/custom/PRJNA512027/metadata/SraRunTable.csv",
    "Chen (GSE213621)": "Human/Patient_Cohorts/pipelines/custom/GSE213621/metadata/SraRunTable.csv",
    "Kozumi (GSE167523)": "Human/Patient_Cohorts/pipelines/custom/GSE167523/metadata/samples.tsv",
    "Hoang (GSE130970)": "Human/Patient_Cohorts/archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv",
    "Suppli (GSE126848)": "Human/Patient_Cohorts/pipelines/custom/GSE126848/metadata/SraRunTable.csv",
    # Mouse
    "FPC/CDAHFD (GSE162876)": "Mouse/Public_Diet_Models/GSE162876/metadata/GSE162876_sample_attributes.csv",
    "LIDPAD (GSE159911)": "Mouse/Public_Diet_Models/GSE159911/metadata/GSE159911_sample_attributes.csv",
    "HFD Short (GSE224069)": "Mouse/Public_Diet_Models/GSE224069/metadata/GSE224069_sample_attributes.csv",
    "HFD Long (GSE274914)": "Mouse/Public_Diet_Models/GSE274914/metadata/GSE274914_sample_attributes.csv",
    "Paquette (GSE156918)": "Mouse/Public_MCD/GSE156918/metadata/samples.tsv",
    "GAN (GSE225616)": "Mouse/Public_Diet_Models/GSE225616/metadata/GSE225616_sample_attributes.csv",
    "In-House (Sanjana Lab)": "Mouse/InHouse_MCD/metadata/samples.tsv",
    "Yue MCD (GSE205974)": "Mouse/Public_MCD/GSE205974/metadata/samples.tsv"
}

# Define Keywords for each category
keywords = {
    "Age": ["age", "year"],
    "Sex": ["sex", "gender"],
    "Labels (NAFL/NASH)": ["diagnosis", "disease", "subtype", "histology", "group", "condition"],
    "Fibrosis": ["fibrosis", "fibrosis_stage", "fibrosis stage", "stage"],
    "NAS": ["nas", "nafld_activity_score", "nas_score", "steatosis", "ballooning", "inflammation", "lobular_inflammation"],
    "Timepoints": ["week", "time", "duration", "timepoint"]
}

print(f"{'Dataset':<30} | {'Found File':<5} | {'Age':<5} | {'Sex':<5} | {'Labels':<8} | {'Fibrosis':<8} | {'NAS':<5} | {'Time':<5}")
print("-" * 100)

for name, rel_path in dataset_paths.items():
    full_path = os.path.join(BASE_DIR, rel_path)
    found = "YES" if os.path.exists(full_path) else "NO"
    
    results = {k: "NO" for k in keywords}
    
    if found == "YES":
        try:
            # Detect separator based on extension
            sep = '\t' if full_path.endswith('.tsv') or full_path.endswith('.txt') else ','
            df = pd.read_csv(full_path, sep=sep, nrows=5)
            cols = [c.lower() for c in df.columns]
            
            for cat, words in keywords.items():
                match = any(any(w in c for w in words) for c in cols)
                if match:
                    results[cat] = "YES"
                    
            # Refine NAS check: needs partial components
            # If NAS is NO, check specifically for components
            if results["NAS"] == "NO":
                 has_steatosis = any("steatosis" in c for c in cols)
                 has_inflam = any("inflammation" in c for c in cols)
                 has_balloon = any("ballooning" in c for c in cols)
                 if has_steatosis or has_inflam or has_balloon:
                     results["NAS"] = "PART"

        except Exception as e:
            found = "ERR"
            print(f"Error reading {name}: {e}")

    print(f"{name:<30} | {found:<10} | {results['Age']:<5} | {results['Sex']:<5} | {results['Labels (NAFL/NASH)']:<8} | {results['Fibrosis']:<8} | {results['NAS']:<5} | {results['Timepoints']:<5}")
