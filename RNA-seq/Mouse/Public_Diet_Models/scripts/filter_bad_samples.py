import pandas as pd
import os

BAD_SAMPLES = ["SRR30252151", "SRR30252152", "SRR12883257", "SRR30252154"]
METADATA_FILE = "metadata/samples.tsv"

if os.path.exists(METADATA_FILE):
    df = pd.read_csv(METADATA_FILE, sep="\t")
    initial_count = len(df)
    
    # Filter
    df_clean = df[~df['sample_id'].isin(BAD_SAMPLES)]
    final_count = len(df_clean)
    
    if initial_count != final_count:
        df_clean.to_csv(METADATA_FILE, sep="\t", index=False)
        print(f"Removed {initial_count - final_count} bad samples. New count: {final_count}")
        print(f"Removed samples: {BAD_SAMPLES}")
    else:
        print("No samples removed. Check IDs?")
else:
    print("Metadata file not found!")
