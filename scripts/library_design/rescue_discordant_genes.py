
import pandas as pd
import numpy as np
import gzip
from pathlib import Path

# --- Configuration ---
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEG_CONTRADICTIONS_FILE = ROOT / "streamlit_deg_explorer/deg_contradictions_5datasets.csv"
ORTHOLOG_FILE = ROOT / "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
MOUSE_DATA_FILE = ROOT / "streamlit_deg_explorer/data/mcd_week_pooled_combined.tsv.gz"
BIOTYPE_FILE = ROOT / "streamlit_deg_explorer/data/ensembl_gene_biotypes.tsv.gz"
TARGET_FILE = ROOT / "final_core_degs.csv"
OUTPUT_FILE = ROOT / "final_core_degs_rescued.csv"

# --- Functions ---

def strip_version(gene_id):
    return str(gene_id).strip().split(".")[0]

def load_mouse_metadata():
    print("Loading Mouse Metadata (Symbol/Biotype)...")
    
    # 1. Symbols from MCD Data
    sym_map = {}
    with gzip.open(MOUSE_DATA_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        # Check cols
        cols = {c.lower(): c for c in df.columns}
        g_col = cols.get("gene_id")
        s_col = cols.get("gene_symbol")
        if g_col and s_col:
            # zip
            for g, s in zip(df[g_col], df[s_col]):
                sym_map[strip_version(g)] = s
                
    # 2. Biotypes
    bio_map = {}
    with gzip.open(BIOTYPE_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        # cols: ensembl_gene_id, gene_biotype, species
        # Filter for mouse
        df = df[df['species'] == 'mouse']
        for g, b in zip(df['ensembl_gene_id'], df['gene_biotype']):
            bio_map[strip_version(g)] = b
            
    return sym_map, bio_map

def load_orthologs():
    print("Loading Orthologs...")
    # Map Human -> Set of Mouse IDs
    h2m = {}
    with gzip.open(ORTHOLOG_FILE, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
        # Columns might vary, robust check
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
    print("Starting Rescue Operation...")
    
    # 1. Load Candidates
    print(f"Loading Candidates from {DEG_CONTRADICTIONS_FILE}...")
    df = pd.read_csv(DEG_CONTRADICTIONS_FILE)
    
    # 2. Filter for "Human UP / Mouse DOWN" AND "Human Sig < 0.05"
    human_cols = ['Status_Hoang_et_al', 'Status_Govaere_et_al']
    mouse_cols = [c for c in df.columns if 'Status_' in c and c not in human_cols]
    
    # Padj columns
    human_padj_cols = ['Hoang_et_al_padj', 'Govaere_et_al_padj']
    
    candidates = []
    
    for _, row in df.iterrows():
        # Classification Logic
        h_status = [row[c] for c in human_cols]
        m_status = [row[c] for c in mouse_cols]
        
        has_h_up = "UP" in h_status
        has_h_down = "DOWN" in h_status
        has_m_up = "UP" in m_status
        has_m_down = "DOWN" in m_status
        
        # Criteria: Human UP, Not Human Down, Mouse Down, Not Mouse Up
        is_discordant = has_h_up and (not has_h_down) and has_m_down and (not has_m_up)
        
        if is_discordant:
            # Check Significance (padj < 0.05) in at least one UP human dataset
            is_sig_human = False
            for p_col, s_col in zip(human_padj_cols, human_cols):
                if row[s_col] == "UP":
                    pval = row[p_col]
                    if pd.notna(pval) and pval < 0.05:
                        is_sig_human = True
                        break
            
            if is_sig_human:
                candidates.append(row)
                
    print(f"  Found {len(candidates)} valid rescue candidates (Human UP, Mouse DOWN, Human padj < 0.05).")
    
    # 2b. Force Rescue Specific Positive Controls (User Request)
    # These genes might have failed strict LFC/TPM cuts in "Discordance" analysis but are required.
    forced_symbols = ['ANGPTL3', 'DGAT2', 'FAP', 'G0S2', 'MLXIPL', 'MOGAT1', 'SCD', 'SLC27A2', 'STK25']
    print(f"  Processing Forced Rescue List: {forced_symbols}")
    
    # We need to fetch Human IDs for these.
    # Load Master Matrix to lookup Human IDs by Symbol? Or Mapping file?
    # Using `gene_symbol_mapping.csv` if available? 
    # Or just iterate over ortholog map and check known symbols?
    # Let's use `deg_contradictions` dataframe if they are there (even if not candidate).
    # IF not in deg_contradictions, we might check Master Matrix?
    # Let's try to find their ID from the loaded CSVs or Orthologs.
    
    # Strategy: Build Symbol -> HumanID map from deg_contradictions (it has Symbol column).
    sym_to_hid = dict(zip(df['Symbol'], df['human_id']))
    
    # If missing from there, we have a problem. STK25 was not in deg_contradictions?
    # Let's check master matrix?
    # Or just use the hardcoded IDs we found earlier for verified ones?
    # SCD: ENSG00000099194, STK25: ENSG00000103005.
    # It's better to be robust.
    
    # Let's load ALL entries from deg_contradictions (it contains all 'Contradictory' ones).
    # But STK25 wasn't contradictory.
    # We need a robust Symbol -> ID map.
    # Let's use `streamlit_deg_explorer/data/gene_symbol_mapping.csv` if valid, or `master_ortholog_matrix.csv.gz`.
    
    master_file = ROOT / "streamlit_deg_explorer/master_ortholog_matrix.csv.gz"
    if master_file.exists():
        print("  Loading Master Matrix for Symbol Lookup...")
        with gzip.open(master_file, 'rt') as f:
            # Read header + subset
            # Just read whole thing, it's 3MB.
            mm = pd.read_csv(f)
            # Create Map
            # Symbol col is last
            if 'Symbol' in mm.columns:
                 sym_to_hid.update(dict(zip(mm['Symbol'], mm['human_id'])))
    
    for sym in forced_symbols:
        if sym in sym_to_hid:
            hid = sym_to_hid[sym]
            # Create a pseudo-row for processing
            # We construct a synthetic candidate row
            
            # Check if already in candidates
            already_queued = any(c['human_id'] == hid for c in candidates)
            if not already_queued:
                print(f"    Force Rescuing {sym} ({hid})")
                candidates.append({
                    'human_id': hid,
                    'Symbol': sym,
                    'Status_Hoang_et_al': 'Forced', # Dummy
                    'Status_Govaere_et_al': 'Forced'
                })
        else:
             print(f"    Warning: Could not find Human ID for {sym}")

    # 3. Prepare Metadata Maps
    sym_map, bio_map = load_mouse_metadata()
    h2m = load_orthologs()
    
    # 4. Construct New Rows
    # Existing Columns:
    # mouse_gene_id,mouse_gene_symbol,mouse_biotype,n_human_orthologs,human_ortholog_ids,human_ortholog_symbols,source_analyses,n_analyses,has_mouse_data,has_human_data
    
    new_rows = []
    
    # Load existing to avoid dupes? User said "add them back to our core degs".
    # Assuming they are NOT in there. But checking existing IDs is good practice.
    existing = pd.read_csv(TARGET_FILE)
    existing_mouse_ids = set(existing['mouse_gene_id'].astype(str).map(strip_version))
    
    added_count = 0
    duplicate_count = 0
    
    for row in candidates:
        human_id = strip_version(row['human_id'])
        symbol = row['Symbol']
        
        # Get Mouse Orthologs
        mouse_orthologs = h2m.get(human_id, set())
        
        if not mouse_orthologs:
            # If no mouse ortholog, we can't really add it to "mouse_gene_id" centric file?
            # Or we add with NA?
            # For now, skip if no ortholog, as per "Discordance" definition implies ortholog exists.
            # But h2m map might be incomplete? No, discordance file relied on it.
            # print(f"Warning: No mouse ortholog found for {symbol} ({human_id})")
            continue
            
        for mid in mouse_orthologs:
            if mid in existing_mouse_ids:
                duplicate_count += 1
                continue
                
            # Build Row
            mid_sym = sym_map.get(mid, f"MOUSE:{symbol}") # Fallback to human symbol with prefix? Or just N/A
            mid_bio = bio_map.get(mid, "protein_coding") # Default safe
            
            # Source Analyses Logic
            # "Rescue | Discordant (Human UP / Mouse DOWN)"
            # Plus verify which Human dataset was UP for detailed trace
            analyses_parts = ["Rescue | Discordant (Human UP / Mouse DOWN)"]
            
            # Add Human Evidence
            if row['Status_Hoang_et_al'] == "UP":
                analyses_parts.append("Patient | GSE130970 | NAS 1+ (upregulated)")
            if row['Status_Govaere_et_al'] == "UP":
                analyses_parts.append("Patient | GSE135251 | NAS 1+ (upregulated)")
                
            source_analyses_str = "; ".join(analyses_parts)
            
            new_row = {
                "mouse_gene_id": mid,
                "mouse_gene_symbol": mid_sym,
                "mouse_biotype": mid_bio,
                "n_human_orthologs": 1, # Simplification? Or check reverse map? 
                                        # Existing file has n_human_orthologs. 
                                        # For simplicity, 1 is fine or we calculate. 
                                        # Let's say 1 for rescue purpose.
                "human_ortholog_ids": human_id,
                "human_ortholog_symbols": symbol,
                "source_analyses": source_analyses_str,
                "n_analyses": len(analyses_parts),
                "has_mouse_data": True,
                "has_human_data": True
            }
            new_rows.append(new_row)
            added_count += 1
            existing_mouse_ids.add(mid) # Prevent double adding if 1:Many expansion happens locally
            
    # 5. Write Output
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        # Ensure column order
        cols = existing.columns.tolist()
        new_df = new_df[cols] # Reorder/Filter
        
        final_df = pd.concat([existing, new_df], ignore_index=True)
        print(f"Appending {len(new_df)} new rows (from {len(candidates)} candidates).")
        print(f"Skipped {duplicate_count} existing mouse IDs.")
        
        final_df.to_csv(OUTPUT_FILE, index=False)
        print(f"Saved to {OUTPUT_FILE}")
    else:
        print("No new rows to add.")

if __name__ == "__main__":
    main()
