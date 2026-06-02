import pandas as pd

# Path
csv_path = 'positive_control.csv'
df = pd.read_csv(csv_path)

# Rename Column
df.rename(columns={'Direction_upon_KD': 'Steatosis_Change_upon_KD'}, inplace=True)

# Define Updates (Gene -> Direction)
updates = {
    'PNPLA3': 'Decrease',   # Therapeutic KD of risk variant / I148M
    'HSD17B13': 'Decrease', # Protective loss-of-function
    'FAP': 'Decrease',      # Anti-fibrotic target
    'IRS1': 'Increase',     # Insulin resistance -> Steatosis
    'PPARG': 'Decrease',    # Hepatic PPARG is lipogenic
    'PPARD': 'Increase',    # Loss of oxidation -> Steatosis
    'LGALS3': 'Decrease',   # Inhibition treats fibrosis/NASH
    'G0S2': 'Decrease',     # Inhibition of ATGL inhibitor -> More lipolysis
    'GLP1R': 'Increase',    # Loss of therapeutic signaling
    'GIPR': 'Increase',
    'GCGR': 'Increase',
    'THRB': 'Increase',
    'FGFR1': 'Increase',
    'FGFR2': 'Increase',
    'FGFR3': 'Increase',
    'FGFR4': 'Increase',
    'FGF21': 'Increase'
}

# Apply Updates
for gene, method in updates.items():
    if gene in df['Gene symbol'].values:
        df.loc[df['Gene symbol'] == gene, 'Steatosis_Change_upon_KD'] = method
        # Also update 'Expected effect' if empty or generic
        current_effect = df.loc[df['Gene symbol'] == gene, 'Expected effect of knockdown on steatosis'].iloc[0]
        if 'Inferred' in str(current_effect):
             df.loc[df['Gene symbol'] == gene, 'Expected effect of knockdown on steatosis'] = f"Loss of target ({method}s steatosis)"

# Save
df.to_csv('positive_control_verified.csv', index=False)
print("Updated CSV saved to positive_control_verified.csv")
print(df[df['Gene symbol'].isin(updates.keys())][['Gene symbol', 'Steatosis_Change_upon_KD']].head(10))
