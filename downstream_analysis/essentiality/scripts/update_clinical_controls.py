import pandas as pd
import numpy as np

# File path
csv_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/positive_control.csv'

# Load existing data
df = pd.read_csv(csv_path)

# Clinical Data from Screenshot
# Structure: Drug, Gene Target (Raw), Company, Phase, Status, Mechanism
clinical_data = [
    {'Drug': 'Tirzepatide', 'Target': 'GLP1R, GIPR', 'Company': 'Eli Lilly', 'Phase': '3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Semaglutide', 'Target': 'GLP1R', 'Company': 'Novo Nordisk', 'Phase': '3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Survodutide', 'Target': 'GLP1R, GCGR', 'Company': 'Boehringer', 'Phase': '3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Efruxifermin', 'Target': 'FGFR1c/2c/3c', 'Company': 'Akero', 'Phase': '3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Pegozafermin', 'Target': 'FGFR1c/2c/3c', 'Company': '89bio', 'Phase': '', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Efimosfermin', 'Target': 'FGFR', 'Company': 'Boston Pharma', 'Phase': '2b', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Resmetirom (Rezdiffra)', 'Target': 'THRB', 'Company': 'Madrigal', 'Phase': '3', 'Status': 'Approved for MASH', 'Mechanism': ''},
    {'Drug': 'Denifanstat', 'Target': 'FASN', 'Company': 'Sagimet', 'Phase': '3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'ION224', 'Target': 'DGAT2', 'Company': 'Ionis/Roche', 'Phase': '2b', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Lanifibranor', 'Target': 'PPARA/D/G', 'Company': 'Inventiva', 'Phase': '3', 'Status': '', 'Mechanism': 'activate PPARs to regulate inflammation...'},
    {'Drug': 'Belapectin', 'Target': 'LGALS3', 'Company': '', 'Phase': '2b/3', 'Status': '', 'Mechanism': ''},
    {'Drug': 'ALN-HSD', 'Target': 'HSD17B13', 'Company': 'Alnylam/Regeneron', 'Phase': '2', 'Status': '', 'Mechanism': 'siRNA'},
    {'Drug': 'AZD7503', 'Target': 'HSD17B13', 'Company': 'AstraZeneca', 'Phase': '1', 'Status': '', 'Mechanism': 'Antisense Oligo (ASO)'},
    {'Drug': 'ALN-PNP', 'Target': 'PNPLA3', 'Company': 'Alnylam/Regeneron', 'Phase': '1', 'Status': '', 'Mechanism': 'siRNA'},
    {'Drug': 'Efruxifermin', 'Target': 'FGF21 Analog', 'Company': '89bio', 'Phase': '2b', 'Status': '', 'Mechanism': 'F4 fibrosis cohort'},
    {'Drug': 'AZD2389', 'Target': 'FAP', 'Company': 'AstraZeneca', 'Phase': '1', 'Status': '', 'Mechanism': 'anti-fibrotic'},
    {'Drug': 'ALG-055009', 'Target': 'THRB', 'Company': 'Aligos', 'Phase': '2', 'Status': '', 'Mechanism': ''},
    {'Drug': 'Various', 'Target': 'IRS1', 'Company': 'Insitro', 'Phase': 'Various', 'Status': '', 'Mechanism': ''}
]

# Create Clinical DataFrame
clin_df = pd.DataFrame(clinical_data)

# Normalize Gene Targets for Mapping
# We need to map 'Target' from clinical to 'Gene symbol' in main df.
# Some targets are comma separated or have slashes.
# We will explode them.

# Expand rows with multiple targets separated by comma
clin_df['Clean_Target'] = clin_df['Target'].str.split(',')
clin_df = clin_df.explode('Clean_Target')
clin_df['Clean_Target'] = clin_df['Clean_Target'].str.strip()

# Handle standard gene mappings
# FGFR1c/2c/3c -> Likely FGFR1, FGFR2, FGFR3
# PPARA/D/G -> PPARA, PPARD, PPARG
# FGF21 Analog -> FGF21

mapping_replacements = {
    'GLP1R': ['GLP1R'],
    'GIPR': ['GIPR'],
    'GCGR': ['GCGR'],
    'FGFR1c/2c/3c': ['FGFR1', 'FGFR2', 'FGFR3'],
    'FGFR': ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'], # Generic?
    'FGF21 Analog': ['FGF21'],
    'PPARA/D/G': ['PPARA', 'PPARD', 'PPARG'],
}

expanded_rows = []
for _, row in clin_df.iterrows():
    t = row['Clean_Target']
    if t in mapping_replacements:
        for new_t in mapping_replacements[t]:
            new_row = row.copy()
            new_row['Clean_Target'] = new_t
            expanded_rows.append(new_row)
    else:
        # Handle 'LGALS3 (Galectin' cleaning
        if 'LGALS3' in t: new_row = row.copy(); new_row['Clean_Target'] = 'LGALS3'; expanded_rows.append(new_row)
        else: expanded_rows.append(row)

clin_df_expanded = pd.DataFrame(expanded_rows)

# Aggregate Clinical Info by Gene
# A gene might have multiple drugs. We concat them.
def agg_clin(x):
    return '; '.join(x.dropna().unique())

clin_agg = clin_df_expanded.groupby('Clean_Target').agg({
    'Drug': agg_clin,
    'Company': agg_clin,
    'Phase': agg_clin,
    'Status': agg_clin,
    'Mechanism': agg_clin
}).reset_index()

# Rename columns for merge
clin_agg.columns = ['Gene symbol', 'Clinical_Drug', 'Clinical_Company', 'Clinical_Phase', 'Clinical_Status', 'Clinical_Mechanism']

# Merge with Main DF
merged_df = pd.merge(df, clin_agg, on='Gene symbol', how='outer')

# Fill NaNs
for c in ['Clinical_Drug', 'Clinical_Company', 'Clinical_Phase', 'Clinical_Status', 'Clinical_Mechanism']:
    merged_df[c] = merged_df[c].fillna('')

# Logic for New Genes (Direction_upon_KD inference or default)
# If a gene was added (only in clinical), 'Direction_upon_KD' will be NaN.
# We can try to infer or fill with 'Unknown'.
# User mainly asked to incorporate.

# Inference Logic based on Drug Mechanism if possible?
# GLP1R (Agonist) -> KD typically Increase Steatosis
# THRB (Agonist) -> KD typically Increase Steatosis
# HSD17B13 (siRNA/ASO) -> KD typically Decrease Steatosis
# PNPLA3 (siRNA) -> KD typically Decrease Steatosis (Existing is Increase, let's respect existing if present)

def infer_direction(row):
    if pd.notna(row['Direction_upon_KD']):
        return row['Direction_upon_KD']
    
    gene = row['Gene symbol']
    mech = row['Clinical_Mechanism']
    drug = row['Clinical_Drug']
    
    # HSD17B13
    if gene == 'HSD17B13': return 'Decrease' # siRNA/ASO treats
    # THRB
    if gene == 'THRB': return 'Increase' # Agonist treats
    # GLP1R
    if gene == 'GLP1R': return 'Increase' # Agonist treats
    # GCGR
    if gene == 'GCGR': return 'Increase' # Agonist treats
    # GIPR
    if gene == 'GIPR': return 'Increase' # Agonist treats
    # FGFRs
    if 'FGFR' in gene: return 'Increase' # Agonist (FGF21 analogs etc) treats
    # LGALS3
    if gene == 'LGALS3': return 'Decrease' # Inhibitor (Belapectin) treats -> KD Decrease
    # FAP
    if gene == 'FAP': return 'Decrease' # Anti-fibrotic (Inhibitor?) -> KD Decrease
    # IRS1
    if gene == 'IRS1': return 'Increase' # Loss usually bad (Insulin Resistance)
    
    return 'TBD' 

merged_df['Direction_upon_KD'] = merged_df.apply(infer_direction, axis=1)

# Sort: Put Clinical matches on top? Or keep original order?
# Let's keep original order (known from First Screenshot) then updated/new ones.
# We can check if 'Pathway / role' is empty for new ones and fill it with 'Target of [Drug]'?
def fill_pathway(row):
    if pd.isna(row['Pathway / role']) or row['Pathway / role'] == '':
        return f"Clinical Target ({row['Clinical_Drug']})"
    return row['Pathway / role']

merged_df['Pathway / role'] = merged_df.apply(fill_pathway, axis=1)
merged_df['Expected effect of knockdown on steatosis'] = merged_df['Expected effect of knockdown on steatosis'].fillna('Inferred from clinical trial mechanism')

# Save
merged_df.to_csv(csv_path, index=False)
print(f"Updated {csv_path} with clinical data.")
print("New genes added:", len(merged_df) - len(df))
print(merged_df[merged_df['Clinical_Drug'] != ''][['Gene symbol', 'Clinical_Drug', 'Clinical_Phase']].head(10))
