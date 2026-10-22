"""Fig 4A label check: printed census numbers vs recomputed atlas census (r4a_atlas_library_counts.tsv, donor roster,
program registry). Also text colours and shapes. Usage: python g4a_check.py <pdf> <tag>"""
import sys, re, subprocess
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
A = pd.read_csv(f"{out}/r4a_atlas_library_counts.tsv", sep="\t")
roster = pd.read_csv(f"{root}/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv", sep="\t")
reg = pd.read_csv(f"{root}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv", sep="\t")
exp = {"datasets": A.dataset.nunique(), "donors_descriptive": len(roster) + 2, "cells_M": round(A.n.sum() / 1e6, 2),
       "analyzed_donors": len(roster), "programs": len(reg), "lineages": reg.cell_type.nunique()}
txt = subprocess.run(["pdftotext", "-layout", pdfp, "-"], capture_output=True, text=True).stdout
flat = " ".join(txt.split())
print("printed text:", flat)
pat = {"datasets": r"(\d+)\s*datasets", "donors_descriptive": r"(\d+)\s*donors", "cells_M": r"([\d.]+)M\s*cells",
       "analyzed_donors": r"(\d+)\s*analyzed donors", "programs": r"(\d+)\s*programs", "lineages": r"(\d+)\s*lineages"}
res = []
for k, p in pat.items():
    m = re.search(p, flat)
    v = float(m.group(1)) if m else None
    res.append(dict(tag=tag, item=k, printed=v, expected=exp[k], off=(v is None) or abs(v - exp[k]) > 1e-9))
R = pd.DataFrame(res); print(R.to_string())
pdf, pg, ch, w, r, l, c = load(pdfp)
print("text colours:", ch.fill.value_counts().to_dict(), "| sizes", ch["size"].round(2).unique().tolist())
print("shapes: rects", len(r), "curves", len(c), "lines", len(l), "| fills", pd.concat([r.fill, c.fill]).value_counts().head(8).to_dict())
R.to_csv(f"{out}/g4a_{tag}.tsv", sep="\t", index=False)
