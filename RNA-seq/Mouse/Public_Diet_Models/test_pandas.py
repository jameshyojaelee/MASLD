import pandas as pd
try:
    idx = pd.Index(["A", "B", "C"])
    meta = ["A"]
    print(f"Index: {idx}")
    print(f"List: {meta}")
    res = idx + meta
    print(f"Result: {res}")
except Exception as e:
    print(f"Error: {e}")
