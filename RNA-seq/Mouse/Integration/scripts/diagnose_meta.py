import pandas as pd
import os

BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
OTHER_META = os.path.join(BASE_DIR, "Mouse/Public_Diet_Models/metadata/samples.tsv")

df = pd.read_csv(OTHER_META, sep="\t")
print("Unique Conditions in Public Diet Meta:")
print(df['condition'].unique())

print("\nChecking GSE224069 specifically:")
print(df[df['dataset'] == 'GSE224069']['condition'].unique())
