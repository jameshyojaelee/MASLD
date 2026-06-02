import pandas as pd

# Load the data
file_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'
df = pd.read_csv(file_path)

# Filter for "not screened in liang" genes
# These are likely currently classified as something else, but I need to check the request.
# The user asked: 'what about in "not screened in liang" genes? how many of those are between -1 and -0.5'
# In previous context (which I don't fully have but can infer), "not screened in liang" might be genes that have a DepMap score but NO Liang score, or genes explicitly marked as such.
# Looking at the file content:
# "Essential_Liang_Only" implies screened in Liang.
# "Essential_DepMap_Only" implies essential in DepMap but NOT in Liang (either not screened or not essential).
# "Essential_Both" implies screened and essential in both.
# "Non_Essential" implies not essential in either (or not screened in one and not essential in other).

# However, the user specifically used the phrase "not screened in liang".
# This usually implies genes present in DepMap (and thus having a DepMap score) but MISSING from the Liang screen data.
# In the `Unified_Essentiality_Combined.csv` (which I viewed earlier), there were columns like `Liang_Pct_Essential` and `DepMap_Pct_Essential`.
# If `Liang_Pct_Essential` is missing (NaN), it might be "not screened".
# But `Essentiality_Quadrants_All.csv` has `Essentiality_Quadrant`.
# Let's assume the user means genes that are essential in DepMap but NOT capable of being assessed in Liang because they weren't screened.
# OR, more simply, genes where `Mean_Essentiality_Score` (which typically comes from Liang) is... wait.
# If they are "not screened in Liang", they wouldn't have a Liang Essentiality Score.
# BUT the User calculates "Mean_Essentiality_Score".
# Let's check `Unified_Essentiality_Combined.csv` again to understand the source of `Mean_Essentiality_Score`.
# If `Mean_Essentiality_Score` is the Liang score, then "not screened in Liang" genes essentially CANNOT have a score between -1 and -0.5 unless the score comes from somewhere else.

# Let's assume the user might be referring to `DepMap_Cas9_Score` for these genes OR `Mean_Essentiality_Score` if it's a composite.
# Actually, looking at `Unified_Essentiality_Combined.csv`:
# Col 2: `Mean_Essentiality_Score`. Col 7: `DepMap_Cas9_Score`.
# If a gene is not screened in Liang, it likely has NaN or 0 for Liang data.
# The user asks for "not screened in liang" genes with "essentiality score between -1 and -0.5".
# If they are not screened, they don't have a Liang score.
# Perhaps the "essentiality score" refers to the DepMap score for these genes?
# OR, perhaps `Mean_Essentiality_Score` is a combined metric? No, usually it's Liang.
# Let's verify if there are genes with `Essentiality_Quadrant` = 'Not_Screened' or similar?
# I saw: Non_Essential, Essential_Both, Essential_DepMap_Only, Essential_Liang_Only.
# Maybe "Essential_DepMap_Only" encompasses "Not Screened in Liang"?
# Or maybe there is a separate list.

# Let's use Python to check for a specific category or patterns in `Unified_Essentiality_Combined.csv` which has more columns.
# I'll check for rows where Liang data is missing but DepMap data exists.

df_unified = pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/results/Unified_Essentiality_Combined.csv')
# Check columns
# Assuming 'Mean_Essentiality_Score' is Liang.
# Assuming 'DepMap_Cas9_Score' is DepMap.

# Let's check if there are genes with DepMap score but NO Mean_Essentiality_Score (or 0/NaN).
# Effectively finding "Not Screened in Liang".
# And if the user asks for "essentiality score", and they are not screened in Liang, they must mean the DepMap score.
# SO: Count genes where (Not Screened in Liang) AND (-1 <= DepMap_Cas9_Score <= -0.5).

# Let's write a script to explore this hypothesis.
print("Columns:", df_unified.columns.tolist())
print("Sample data:\n", df_unified.head())

# Check for missing Liang scores
missing_liang = df_unified[df_unified['Mean_Essentiality_Score'].isna()]
print(f"Genes with missing Liang score: {len(missing_liang)}")

# Check for 0 Liang scores that might indicate not screened
zero_liang = df_unified[df_unified['Mean_Essentiality_Score'] == 0]
print(f"Genes with 0 Liang score: {len(zero_liang)}")

