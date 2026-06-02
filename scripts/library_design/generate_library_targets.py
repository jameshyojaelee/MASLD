
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Definitions
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DATA_DIR = ROOT / "streamlit_deg_explorer" / "data"

INHOUSE_MCD_FILES = {
    "MCD (in-house) | mouse | MCD Week pooled (combined)": "mcd_week_pooled_combined.tsv.gz",
}

EXTERNAL_MCD_FILES = {
    "MCD (external) | mouse | GSE156918 (external MCD)": "other_mcd_gse156918.tsv.gz",
    "MCD (external) | mouse | GSE205974 (external MCD)": "other_mcd_gse205974.tsv.gz",
}

# Note: Using NAS 1+ vs 0 results as the "NAS high/1+" set as per recent verification tasks
PATIENT_FILES = {
    "GSE130970": {
        "na_label": "Patient | GSE130970 | NAS 1+ (upregulated)",
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/analysis/differential_expression/current/nas_threshold_sensitivity/cumulative_nas/GSE130970/nas_1_vs_0/results.csv",
        "fib_label": "Patient | GSE130970 | Fibrosis (upregulated)",
        "fibrosis": "gse130970_fibrosis.csv.gz",
    },
    "GSE135251": {
        "na_label": "Patient | GSE135251 | NAS 1+ (upregulated)",
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/analysis/differential_expression/current/nas_threshold_sensitivity/cumulative_nas/GSE135251/nas_1_vs_0/results.csv",
        "fib_label": "Patient | GSE135251 | Fibrosis (upregulated)",
        "fibrosis": "gse135251_fibrosis.csv.gz",
    },
}

GWAS_FILE = "Closest_genes.csv"
GWAS_LABEL = "GWAS | human | Closest genes"

ORTHOLOG_FILENAME = "mouse_human_orthologs.tsv.gz"
ORTHOLOG_PATH = DATA_DIR / ORTHOLOG_FILENAME

TPM_FILES = {
    "GSE130970": DATA_DIR / "gse130970_nas_high.csv.gz",
    "GSE135251": DATA_DIR / "gse135251_nas_high.csv.gz"
}

GENE_SYMBOL_MAP = DATA_DIR / "gene_symbol_mapping.csv"
GENE_BIOTYPE_MAP = DATA_DIR / "ensembl_gene_biotypes.tsv.gz"

def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]

def load_df(path, tpm_path=None):
    sep = "\t" if str(path).endswith(".tsv.gz") else ","
    df = pd.read_csv(path, sep=sep)
    cols = {c.lower(): c for c in df.columns}
    rename = {}
    if "gene" in cols: rename[cols["gene"]] = "gene_id"
    elif "gene_id" in cols: rename[cols["gene_id"]] = "gene_id"
    
    if "gene_symbol" in cols: rename[cols["gene_symbol"]] = "gene_symbol"
    elif "symbol" in cols: rename[cols["symbol"]] = "gene_symbol"
    
    if "log2foldchange" in cols: rename[cols["log2foldchange"]] = "log2FoldChange"
    if "padj" in cols: rename[cols["padj"]] = "padj"
    
    if "tpm_mean" in cols: rename[cols["tpm_mean"]] = "tpm_mean"
    elif "tpm" in cols: rename[cols["tpm"]] = "tpm_mean"
    
    # Biotype if present
    if "biotype" in cols: rename[cols["biotype"]] = "biotype"
    elif "gene_biotype" in cols: rename[cols["gene_biotype"]] = "biotype"

    df = df.rename(columns=rename)
    
    # Merge TPM if missing
    if "tpm_mean" not in df.columns and tpm_path and os.path.exists(tpm_path):
        tpm_sep = "\t" if str(tpm_path).endswith(".tsv.gz") else ","
        tpm_df = pd.read_csv(tpm_path, sep=tpm_sep)
        t_cols = {c.lower(): c for c in tpm_df.columns}
        t_gene = t_cols.get("gene_id") or t_cols.get("gene")
        t_val = t_cols.get("tpm_mean") or t_cols.get("tpm")
        
        if t_gene and t_val:
            # Map based on stripped version
            tmap = dict(zip(tpm_df[t_gene].astype(str).map(strip_version), tpm_df[t_val]))
            df["tpm_mean"] = df["gene_id"].astype(str).map(strip_version).map(tmap)
            
    return df

def load_ortholog_map(path: Path):
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    mouse_col = cols.get("mouse_ensembl_gene_id") or cols.get("ensembl_gene_id") or cols.get("mouse_gene_id")
    human_col = cols.get("human_ensembl_gene_id") or cols.get("hsapiens_homolog_ensembl_gene")
    
    df["mouse_ensembl_gene_id"] = df[mouse_col].map(strip_version)
    df["human_ensembl_gene_id"] = df[human_col].map(strip_version)
    
    # Human -> Mouse Mapping (One-to-Many support)
    h2m = {}
    # Mouse -> Human Mapping (For checking existence)
    m2h = {}
    
    for row in df.itertuples():
        m_id = row.mouse_ensembl_gene_id
        h_id = row.human_ensembl_gene_id
        h2m.setdefault(h_id, set()).add(m_id)
        m2h.setdefault(m_id, set()).add(h_id)
        
    return h2m, m2h

def load_symbol_map(path: Path):
    df = pd.read_csv(path)
    mapping = dict(zip(df["gene_id"].astype(str).map(strip_version), df["gene_symbol"]))
    sym2id = {}
    for k, v in mapping.items():
        if pd.notna(v):
            sym2id[v] = k
    return mapping, sym2id

def load_biotype_map(path: Path):
    df = pd.read_csv(path, sep="\t")
    cols = df.columns
    id_col = cols[0]
    bio_col = cols[1]
    mapping = dict(zip(df[id_col].astype(str).map(strip_version), df[bio_col]))
    return mapping

def main():
    print("Loading resources... (Hybrid Orthology Mode)")
    
    # --- CONFIGURATION ---
    CONFIG = {
        "padj_cut": 0.1,
        "lfc_cut": 0.8,
        "tpm_cut_default": 1.0, # Protein coding, miRNA
        "tpm_cut_lnc": 0.5,     # lncRNA
        "allowed_biotypes": {
            "protein_coding", 
            "lncRNA", 
            "miRNA"
        }
    }
    print(f"Configuration: {CONFIG}")
    # ---------------------
    
    # --- MANUAL PATCHES ---
    SYMBOL_ALIASES = {
        "Sepp1": "Selenop"
    }
    # Pten -> PTEN, Scd1 -> SCD
    ORTHOLOGY_PATCH = {
        "ENSMUSG00000013663": "ENSG00000171862",
        "ENSMUSG00000037071": "ENSG00000099194"
    }
    # ----------------------

    h2m, m2h = load_ortholog_map(ORTHOLOG_PATH)
    
    # Apply Orthology Patch
    for m_id, h_id in ORTHOLOGY_PATCH.items():
        if m_id not in m2h:
            m2h[m_id] = set()
        m2h[m_id].add(h_id)
        if h_id not in h2m:
            h2m[h_id] = set()
        h2m[h_id].add(m_id)
        print(f"Patched Orthology for {m_id} -> {h_id}")
        
    id2sym, sym2id = load_symbol_map(GENE_SYMBOL_MAP)
    id2bio = load_biotype_map(GENE_BIOTYPE_MAP)
    
    # --- BUILD MOUSE SYMBOL MAP ---
    print("Building Mouse ID->Symbol map from input files...")
    mouse_id2sym = {}
    
    # Helper to scan a file for symbols
    def scan_for_symbols(fpath):
        try:
            temp_df = load_df(DATA_DIR / fpath)
            if "gene_id" in temp_df.columns and "gene_symbol" in temp_df.columns:
                # Create map
                clean_ids = temp_df["gene_id"].astype(str).map(strip_version)
                symbols = temp_df["gene_symbol"].astype(str)
                # Update dict
                mouse_id2sym.update(dict(zip(clean_ids, symbols)))
        except Exception as e:
            print(f"Warning: Could not scan {fpath} for symbols: {e}")

    # Scan in-house and external MCD files
    for fname in INHOUSE_MCD_FILES.values():
        scan_for_symbols(fname)
    for fname in EXTERNAL_MCD_FILES.values():
        scan_for_symbols(fname)
        
    print(f"  Recovered {len(mouse_id2sym)} mouse symbols.")
    
    # Create Inverse Map (Symbol -> ID) for list lookup
    mouse_sym2id = {}
    for mid, sym in mouse_id2sym.items():
        if pd.notna(sym):
            mouse_sym2id[str(sym)] = mid

    # Merge into main id2sym
    id2sym.update(mouse_id2sym)
    # -----------------------------
    
    datasets = {}
    for label, fname in INHOUSE_MCD_FILES.items():
        datasets[label] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
    for label, fname in EXTERNAL_MCD_FILES.items():
        datasets[label] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
    for study_id, info in PATIENT_FILES.items():
        datasets[info["na_label"]] = {"df": load_df(info["nas_1plus"], TPM_FILES[study_id]), "species": "human"}
        datasets[info["fib_label"]] = {"df": load_df(DATA_DIR / info["fibrosis"], TPM_FILES[study_id]), "species": "human"}
    
    # GWAS (Human List)
    gwas_df = pd.read_csv(DATA_DIR / GWAS_FILE, header=None, names=["gene_symbol"])
    datasets[GWAS_LABEL] = {"df": gwas_df, "species": "human", "is_list": True}

    # PERTURB Multimodal (Mouse List)
    perturb_path = ROOT / "Perturb-Multimodal.csv"
    if perturb_path.exists():
        perturb_df = pd.read_csv(perturb_path, header=None, names=["gene_symbol"])
        datasets["Perturb | mouse | Multimodal Prior"] = {"df": perturb_df, "species": "mouse", "is_list": True}
    else:
        print(f"Warning: {perturb_path} not found.")
    
    target_pool = {}

    print("Processing datasets...")
    for label, info in datasets.items():
        df = info["df"]
        candidates = set()
        
        if info.get("is_list"):
            # Gene Lists (No LFC/Padj)
            for raw_sym in df["gene_symbol"].dropna().astype(str).str.strip():
                sym = SYMBOL_ALIASES.get(raw_sym, raw_sym) # Handle Aliases
                
                gid = None
                
                if info["species"] == "human":
                    # Human List (GWAS) -> Map using standard Human Symbol Map
                    if sym in sym2id:
                        gid = sym2id[sym]
                else:
                    # Mouse List (Perturb) -> Map using inferred Mouse Symbol Map
                    if sym in mouse_sym2id:
                        gid = mouse_sym2id[sym]
                
                if gid:
                    # Check Biotype
                    bio = id2bio.get(gid, "unknown")
                    if bio in CONFIG["allowed_biotypes"]:
                        candidates.add(gid)
        else:
            # Transcriptomics (Standard)
            mask_base = (df["padj"] < CONFIG["padj_cut"]) & (df["log2FoldChange"] > CONFIG["lfc_cut"])
            df_sig = df[mask_base].copy()
            
            # Map IDs
            df_sig["clean_id"] = df_sig["gene_id"].astype(str).map(strip_version)
            
            # Biotypes
            biotypes = df_sig["clean_id"].map(id2bio).fillna("unknown")
            is_lnc = (biotypes == "lncRNA")
            is_allowed = biotypes.isin(CONFIG["allowed_biotypes"])
            
            # TPMs
            if "tpm_mean" in df_sig.columns:
                tpms = df_sig["tpm_mean"].fillna(0)
            else:
                tpms = pd.Series([9999.0] * len(df_sig), index=df_sig.index)
            
            # Logic: Must be allowed biotype AND meet TPM cutoff
            mask_lnc_pass = is_lnc & (tpms >= CONFIG["tpm_cut_lnc"])
            mask_other_pass = (~is_lnc) & (tpms >= CONFIG["tpm_cut_default"])
            
            # Final combined mask: (Allowed Biotype) AND ( (lncRNA & HighTPM) OR (Other & HighTPM) )
            # Since mask_lnc/other imply the biotype check logic implicitly but let's be safe
            final_mask = is_allowed & (mask_lnc_pass | mask_other_pass)
            
            candidates = set(df_sig.loc[final_mask, "clean_id"])

        # 2. Map & Filter
        for gene_id in candidates:
            mouse_ids = set()
            
            if info["species"] == "mouse":
                # Mouse Source
                bio = id2bio.get(gene_id, "unknown")
                
                # Check biotype allowed (redundant if filtering worked above, but safe)
                if bio not in CONFIG["allowed_biotypes"]:
                    continue

                if bio in ["lncRNA", "miRNA"]:
                    mouse_ids.add(gene_id)
                else:
                    # Enforce Orthology for Protein Coding
                    if gene_id in m2h:
                        mouse_ids.add(gene_id)
                    else:
                        continue
                human_ref = None 
                
            else:
                # Human Source
                if gene_id in h2m:
                    mouse_ids = h2m[gene_id]
                    human_ref = gene_id
                else:
                    continue
            
            # 3. Add to Pool
            for mid in mouse_ids:
                # Check Mouse Biotype (for Human -> Mouse mapping cases)
                # Ensure the mapped mouse gene is also an allowed biotype!
                m_bio = id2bio.get(mid, "unknown")
                if m_bio not in CONFIG["allowed_biotypes"]:
                    continue
                    
                if mid not in target_pool:
                    target_pool[mid] = {
                        "human_orthologs": set(),
                        "analyses": set(),
                        "species_source": set()
                    }
                
                target_pool[mid]["analyses"].add(label)
                target_pool[mid]["species_source"].add(info["species"])
                
                if human_ref:
                    target_pool[mid]["human_orthologs"].add(human_ref)
                
                if info["species"] == "mouse":
                    if mid in m2h:
                        target_pool[mid]["human_orthologs"].update(m2h[mid])

    print(f"Total Mouse Targets Investigated: {len(target_pool)}")
    
    # 4. Generate Final Rows
    final_rows = []
    for mid, data in target_pool.items():
        m_bio = id2bio.get(mid, "unknown")
        m_sym = id2sym.get(mid, "") 
        h_orthologs = sorted(list(data["human_orthologs"]))
        h_symbols = [id2sym.get(h, "") for h in h_orthologs]
        analyses = sorted(list(data["analyses"]))
        
        row = {
            "mouse_gene_id": mid,
            "mouse_gene_symbol": m_sym,
            "mouse_biotype": m_bio,
            "n_human_orthologs": len(h_orthologs),
            "human_ortholog_ids": ";".join(h_orthologs),
            "human_ortholog_symbols": ";".join(h_symbols),
            "source_analyses": "; ".join(analyses),
            "n_analyses": len(analyses),
            "has_mouse_data": "mouse" in data["species_source"],
            "has_human_data": "human" in data["species_source"]
        }
        final_rows.append(row)
        
    df_all = pd.DataFrame(final_rows)
    df_all = df_all.sort_values(["n_analyses", "mouse_gene_symbol"], ascending=[False, True])
    
    outfile_core = ROOT / "final_core_degs.csv"
    df_all.to_csv(outfile_core, index=False)
    print(f"Saved {len(df_all)} total core targets to {outfile_core}")
    
    df_lnc = df_all[df_all["mouse_biotype"] == "lncRNA"]
    outfile_lnc = ROOT / "final_lncrna_degs.csv"
    df_lnc.to_csv(outfile_lnc, index=False)
    print(f"Saved {len(df_lnc)} lncRNA targets to {outfile_lnc}")

if __name__ == "__main__":
    main()
