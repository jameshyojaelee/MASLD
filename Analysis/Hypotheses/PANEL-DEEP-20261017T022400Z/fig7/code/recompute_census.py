"""7A census recount from the current exposure_audit.json files (the panel states a September 1 snapshot;
files modified since then will differ). Counts checkpoint findings by exposure_state."""
import json, glob, collections, sys, os
import pandas as pd
M = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/config/artifacts/models"
dirs = sorted(d for d in os.listdir(M) if os.path.isdir(os.path.join(M, d)))
rows = []
for d in dirs:
    f = os.path.join(M, d, "exposure_audit.json")
    if not os.path.exists(f):
        rows.append({"model_dir": d, "finding": None, "exposure_state": "NO_AUDIT_FILE", "mtime": None}); continue
    j = json.load(open(f)); mt = os.path.getmtime(f)
    if "checkpoint_findings" in j and isinstance(j["checkpoint_findings"], dict):
        for k, v in j["checkpoint_findings"].items():
            st = v.get("exposure_state") if isinstance(v, dict) else None
            rows.append({"model_dir": d, "finding": k, "exposure_state": st, "mtime": mt, "form": "plural"})
    if "checkpoint_finding" in j:
        v = j["checkpoint_finding"]; st = v.get("exposure_state") if isinstance(v, dict) else None
        rows.append({"model_dir": d, "finding": "(single)", "exposure_state": st, "mtime": mt, "form": "single"})
    if "checkpoint_finding" not in j and "checkpoint_findings" not in j:
        rows.append({"model_dir": d, "finding": None, "exposure_state": "NO_FINDING", "mtime": mt, "form": "none"})
t = pd.DataFrame(rows)
t["mtime_iso"] = pd.to_datetime(t.mtime, unit="s").dt.strftime("%Y-%m-%d")
t.to_csv(sys.argv[1] + "/census_recount.tsv", sep="\t", index=False)
print("dirs", len(dirs), "audit files", sum(os.path.exists(os.path.join(M, d, "exposure_audit.json")) for d in dirs))
print(t.groupby("form", dropna=False).size())
print(t.exposure_state.value_counts(dropna=False).to_string())
print("files modified after 2026-09-01:", t[t.mtime_iso > "2026-09-01"].model_dir.nunique())
