
import pandas as pd
import numpy as np
import gzip
from pathlib import Path

# --- Configuration ---
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MASTER_MATRIX_FILE = ROOT / "streamlit_deg_explorer/master_ortholog_matrix.csv.gz"
ORTHOLOG_FILE = ROOT / "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
MOUSE_DATA_FILE = ROOT / "streamlit_deg_explorer/data/mcd_week_pooled_combined.tsv.gz"
BIOTYPE_FILE = ROOT / "streamlit_deg_explorer/data/ensembl_gene_biotypes.tsv.gz"
TARGET_FILE = ROOT / "final_core_degs.csv"
OUTPUT_FILE = ROOT / "final_core_degs_updated.csv"

# Strict Rule Implementation (User Mandate)
# Rules: 
# 1. log2FC > 0.8 AND padj < 0.1
# 2. TPM > 1.0 (PCG) OR TPM > 0.5 (lncRNA)
# 3. Strict Orthology (Must have Mouse ID)
# 4. Include ANY Human UP regardless of Mouse Concordance

LFC_CUT = 0.8
PADJ_CUT = 0.1



def strip_version(gene_id):
    if pd.isna(gene_id): return ""
    return str(gene_id).strip().split(".")[0]

def load_human_biotypes():
    print("Loading Human Biotypes...")
    bio_map = {}
    with gzip.open(BIOTYPE_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        df = df[df['species'] == 'human']
        for g, b in zip(df['ensembl_gene_id'], df['gene_biotype']):
            bio_map[strip_version(g)] = b
    return bio_map

def load_mouse_metadata():
    print("Loading Mouse Metadata (Symbol/Biotype)...")
    
    # 1. Symbols from MCD Data
    sym_map = {}
    with gzip.open(MOUSE_DATA_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        cols = {c.lower(): c for c in df.columns}
        g_col = cols.get("gene_id")
        s_col = cols.get("gene_symbol")
        if g_col and s_col:
            for g, s in zip(df[g_col], df[s_col]):
                sym_map[strip_version(g)] = s
                
    # 2. Biotypes
    bio_map = {}
    with gzip.open(BIOTYPE_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        df = df[df['species'] == 'mouse']
        for g, b in zip(df['ensembl_gene_id'], df['gene_biotype']):
            bio_map[strip_version(g)] = b
            
    return sym_map, bio_map

def load_orthologs():
    print("Loading Orthologs...")
    h2m = {}
    with gzip.open(ORTHOLOG_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        cols = {c.lower(): c for c in df.columns}
        h_col = cols.get('human_ensembl_gene_id') or 'human_ensembl_gene_id'
        m_col = cols.get('mouse_ensembl_gene_id') or 'mouse_ensembl_gene_id'
        
        for h, m in zip(df[h_col], df[m_col]):
             if pd.notna(h) and pd.notna(m):
                 h_clean = strip_version(h)
                 m_clean = strip_version(m)
                 if h_clean not in h2m: h2m[h_clean] = set()
                 h2m[h_clean].add(m_clean)
    return h2m

def main():
    print("Starting Library Update (Strict Rules)...")
    
    # Load Biotypes first for filtering
    human_biotypes = load_human_biotypes()
    
    # 1. Load Master Matrix
    print(f"Loading Master Matrix from {MASTER_MATRIX_FILE}...")
    with gzip.open(MASTER_MATRIX_FILE, 'rt') as f:
        master = pd.read_csv(f)
        
    print(f"  Total Rows in Matrix: {len(master)}")
    
    # 2. Identify Human UP Genes (Strict Rules)
    
    candidates = []
    
    for idx, row in master.iterrows():
        human_id = strip_version(row['human_id'])
        symbol = row.get('Symbol', 'N/A')
        biotype = human_biotypes.get(human_id, "protein_coding") # Default to PCG strictness if unknown
        
        # Determine TPM Filter
        if biotype == "lncRNA":
            tpm_cut = 0.5
        else:
            tpm_cut = 1.0
            
        sources = []
        
        # Helper check
        def check_dataset(prefix, label):
            lfc = row.get(f"{prefix}_lfc", -99)
            padj = row.get(f"{prefix}_padj", 1.0)
            tpm = row.get(f"{prefix}_tpm", 0.0)
            
            if pd.isna(lfc) or pd.isna(padj) or pd.isna(tpm): return False
            return (lfc > LFC_CUT) and (padj < PADJ_CUT) and (tpm > tpm_cut)

        if check_dataset("Hoang_et_al", "Hoang"):
            sources.append("Patient | GSE130970 | NAS 1+ (upregulated)")
            
        if check_dataset("Govaere_et_al", "Govaere"):
            sources.append("Patient | GSE135251 | NAS 1+ (upregulated)")
            
        if sources:
            candidates.append({
                "human_id": human_id,
                "symbol": symbol,
                "sources": sources
            })
            
    print(f"  Found {len(candidates)} Human UP genes (Strict Rules).")
    
    # 3. Load Mouse Metadata & Orthologs
    sym_map, bio_map = load_mouse_metadata()
    h2m = load_orthologs()
    
    # 4. Filter Against Existing Library
    print(f"Loading Existing Library from {TARGET_FILE}...")
    existing = pd.read_csv(TARGET_FILE)
    existing_mouse_ids = set(existing['mouse_gene_id'].astype(str).map(strip_version))
    
    new_rows = []
    skipped_count = 0
    added_count = 0
    
    for cand in candidates:
        hid = cand['human_id']
        sym = cand['symbol']
        sources = cand['sources']
        
        m_ids = h2m.get(hid, set())
        
        if not m_ids:
            # Strict Orthology Rule - Skip if no mouse ID
            continue
            
        for mid in m_ids:
            if mid in existing_mouse_ids:
                skipped_count += 1
                continue
                
            # Create New Row
            mid_sym = sym_map.get(mid, f"MOUSE:{sym}")
            mid_bio = bio_map.get(mid, "protein_coding")
            
            source_str = "; ".join(sources)
            
            new_row = {
                "mouse_gene_id": mid,
                "mouse_gene_symbol": mid_sym,
                "mouse_biotype": mid_bio,
                "n_human_orthologs": 1,
                "human_ortholog_ids": hid,
                "human_ortholog_symbols": sym,
                "source_analyses": source_str,
                "n_analyses": len(sources),
                "has_mouse_data": True,
                "has_human_data": True
            }
            new_rows.append(new_row)
            added_count += 1
            existing_mouse_ids.add(mid)
            
    print(f"  Skipped {skipped_count} existing mouse IDs.")
    print(f"  Adding {added_count} new entries.")
    
    # 6. Write Output
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        cols = existing.columns.tolist()
        for c in cols:
            if c not in new_df.columns:
                new_df[c] = None 
        new_df = new_df[cols]
        
        final_df = pd.concat([existing, new_df], ignore_index=True)
        final_df.to_csv(OUTPUT_FILE, index=False)
        print(f"Saved updated library to {OUTPUT_FILE}")
    else:
        print("No new entries to add.")

if __name__ == "__main__":
    main()
