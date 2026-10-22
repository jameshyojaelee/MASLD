import pandas as pd
pd.set_option('display.width',300); pd.set_option('display.max_columns',40); pd.set_option('display.max_colwidth',80)
d=pd.read_parquet('passport_program_context.parquet')
print(d.shape); print(d.columns.tolist())
print(d.drop(columns=[c for c in d.columns if 'sha' in c or 'node' in c]).to_string())
