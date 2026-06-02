import pandas as pd
import os

# Paths
base_dir = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
pos_ctrl_path = os.path.join(base_dir, 'positive_control.csv')
core_degs_path = os.path.join(base_dir, 'final_core_degs.csv')

# Load Data
print(f"Loading {pos_ctrl_path}...")
df_pos = pd.read_csv(pos_ctrl_path)
print(f"Loading {core_degs_path}...")
df_core = pd.read_csv(core_degs_path)

# Prepare Core DEG Symbols
core_symbols = set()
for list_str in df_core['human_ortholog_symbols'].dropna():
    parts = str(list_str).split(';')
    for p in parts:
        core_symbols.add(p.strip())

print(f"Total unique human orthologs in Core DEGs: {len(core_symbols)}")

# Add Column
df_pos['In_Final_Library'] = df_pos['Gene symbol'].apply(lambda x: x in core_symbols)

# Save
df_pos.to_csv(pos_ctrl_path, index=False)
print(f"Updated {pos_ctrl_path} with 'In_Final_Library' column.")
print(df_pos[['Gene symbol', 'In_Final_Library']].head())
