#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --job-name=resolution_r8
# Tests, then the r7 resolution with the 855 unmatched annotations held as not testable.
set -euo pipefail
root=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
src=$root/scripts/analysis/alt_contig_recount
base=$root/figures/candidates/resource-corrections-20261005-v1
followup=$base/bulk_execution_sources/followup
r7=$base/annotation_resolution_r7
r8=$base/annotation_resolution_r8
export PATH="$bin:$PATH"
cd "$src"
"$bin/python" -B -m unittest test_annotation_correspondence test_verified_fragments test_count_execution_evidence test_native_library_subsample
cd "$root"
"$bin/Rscript" --vanilla scripts/analysis/alt_contig_recount/test_assemble_verified_counts.R
"$bin/Rscript" --vanilla scripts/manuscript/resource_f_five/test_qualified_bulk_inputs.R
"$bin/python" -B "$src/resolve_annotation_correspondence.py" \
  --gtf "$(${bin}/python -c "import json;print(json.load(open('$r7/resolution_summary.json'))['gtf'])")" \
  --output "$r8" --source-cache-manifest "$r7/source_manifest.tsv" \
  --prior-correspondence "$r7/annotation_correspondence.tsv" \
  --ensembl-source-root "$followup" --added-patch-source-root "$followup" \
  --grc-exclusions "$followup/grc_exclusions.bed" --refseq-transcript-metadata "$followup/gencode49_RefSeq.gz" \
  --apply-grc-reference-exclusions --hold-unresolved-as-untestable
"$bin/python" -B - "$r7" "$r8" <<'PY'
import csv, gzip, json, re, sys
from collections import Counter
r7, r8 = sys.argv[1:]
read = lambda root: {r["gene_id"]: r for r in csv.DictReader(open(root+"/annotation_correspondence.tsv"), delimiter="\t")}
old, new = read(r7), read(r8)
assert list(old) == list(new), "annotation set or order differs"
changed = Counter()
for gid, a in old.items():
    b = new[gid]
    # Source files are copied per run, so their recorded directory differs by construction.
    fields = tuple(k for k in a if a[k].replace(r7, r8) != b[k])
    if fields: changed[(a["disposition"], b["disposition"], fields)] += 1
summary = json.load(open(r8+"/resolution_summary.json"))
targets = Counter()
with gzip.open(summary["grouped_counting_annotation"], "rt") as handle:
    for line in handle:
        if not line.startswith("#"): targets[re.search(r'gene_id "([^"]+)"', line).group(1)] += 1
held = {g for g, r in new.items() if r["disposition"] == "held_untestable"}
named = {r["gene_name"]: (r["gene_id"], r["disposition"]) for r in new.values()
         if r["gene_name"] in {"APOB", "COL3A1", "PKLR", "TDO2", "GLUD2"}}
report = {"differences_r7_to_r8": {" | ".join(map(str, k)): v for k, v in changed.items()},
          "dispositions": summary["dispositions"], "n_unresolved_annotations": summary["n_unresolved_annotations"],
          "n_held_untestable_annotations": summary["n_held_untestable_annotations"],
          "counting_targets_in_grouped_annotation": len(targets), "held_targets_in_grouped_annotation": len(held & set(targets)),
          "held_biotypes": dict(Counter(new[g]["gene_biotype"] for g in held).most_common()),
          "named_genes_last_row": named,
          "named_gene_names_held": sorted({new[g]["gene_name"] for g in held} & {"APOB", "COL3A1", "PKLR", "TDO2", "GLUD2"})}
only_relabel = set(changed) == {("unresolved", "held_untestable", ("disposition",))} and sum(changed.values()) == 855
report["status"] = "PASS" if only_relabel and not report["named_gene_names_held"] else "FAIL"
json.dump(report, open(r8+"/review_checks.json", "x"), indent=2)
print(json.dumps(report, indent=2))
assert report["status"] == "PASS"
PY
"$bin/python" -m pip freeze > "$r8/pip_freeze.txt"
