#!/usr/bin/env python3
"""
Classify per-dataset strandedness from the mouse strand-check run_info.json files.
For a stranded library, one of rf/fr keeps ~all reads (~unstranded) and the other
~0; an unstranded library keeps ~half in each. Rule on ratios to unstranded:
  rf>0.70 & fr<0.40  -> reverse-stranded (rf)
  fr>0.70 & rf<0.40  -> forward-stranded (fr)
  else               -> unstranded
Writes strand_map_mouse.tsv (dataset \t strand \t p_unstr \t p_rf \t p_fr).
"""
import glob, json, os, sys

ID = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/isoform_diversity"
SC = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
      "RNA-seq/results/isoform_diversity/mouse/_strandcheck")

def p(ds, mode):
    fp = f"{SC}/{ds}_{mode}/run_info.json"
    try:
        return float(json.load(open(fp))["p_pseudoaligned"])
    except Exception:
        return None

datasets = sorted({os.path.basename(d).rsplit("_", 1)[0]
                   for d in glob.glob(f"{SC}/*_unstranded")})
rows = []
for ds in datasets:
    pu, pr, pf = p(ds, "unstranded"), p(ds, "rf"), p(ds, "fr")
    if not pu:
        rows.append((ds, "unstranded", pu, pr, pf)); continue
    rr, rf_ = (pr or 0) / pu, (pf or 0) / pu
    if rr > 0.70 and rf_ < 0.40:
        s = "rf"
    elif rf_ > 0.70 and rr < 0.40:
        s = "fr"
    else:
        s = "unstranded"
    rows.append((ds, s, pu, pr, pf))

out = f"{ID}/strand_map_mouse.tsv"
with open(out, "w") as fh:
    fh.write("dataset\tstrand\tp_unstr\tp_rf\tp_fr\n")
    for ds, s, pu, pr, pf in rows:
        fh.write(f"{ds}\t{s}\t{pu}\t{pr}\t{pf}\n")
        sys.stderr.write(f"[strand] {ds}: {s}  (unstr={pu} rf={pr} fr={pf})\n")
sys.stderr.write(f"[strand] wrote {out}\n")
