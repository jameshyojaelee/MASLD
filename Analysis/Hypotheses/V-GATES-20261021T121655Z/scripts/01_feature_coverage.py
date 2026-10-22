"""WP0a: coverage of the 117 frozen programs from feature or panel files only.

No expression values are read. Presence is a gene on the assay's feature axis;
the detection filter (>=1% of units) is applied later, inside each dataset's
own gate, so presence-only weight is an upper bound on retained weight.

Usage: python -I 01_feature_coverage.py <run_dir>
"""
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import h5py

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
HOTSPOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
REGISTRY_SHA = "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
GENCODE = ROOT / "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"
CONFIRMED = "gencode_v49_unique_symbol_confirmed"
CONFIRMATORY = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "Stromal ECM (IGFBP7)",
    "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular injury (BICC1)",
}
# Spatial observability rule (spatial_resource_lib.py) and single-cell replication rule (H90 G3).
OBS_MIN_GENES, OBS_MIN_WEIGHT = 8, 0.20
REP_MIN_GENES, REP_MIN_WEIGHT = 8, 0.80
CONTROL_PREFIXES = ("NegPrb", "Negative", "SystemControl", "FalseCode", "falsecode", "NegControl")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_programs():
    if sha256(HOTSPOT / "program_registry_v2.tsv") != REGISTRY_SHA:
        raise SystemExit("registry sha256 mismatch")
    with open(HOTSPOT / "program_registry_v2.tsv") as fh:
        registry = list(csv.DictReader(fh, delimiter="\t"))
    if len(registry) != 117:
        raise SystemExit("expected 117 programs")
    members = defaultdict(list)
    with open(HOTSPOT / "program_membership_v2.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["mapped_symbol_status"] == CONFIRMED:
                members[row["program_uid"]].append((row["mapped_symbol"], float(row["original_l1_weight"])))
    return registry, members


def load_gencode():
    sym_to_ens, biotype = defaultdict(set), {}
    with gzip.open(GENCODE, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym_to_ens[row["gene_name"]].add(row["ensembl_base"])
            biotype[row["gene_name"]] = row["gene_biotype"]
    return sym_to_ens, biotype


def axis_10x_features(path):
    ens, sym = set(), set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[2] == "Gene Expression":
                ens.add(f[0].split(".")[0])
                sym.add(f[1])
    return ens, sym


def axis_10x_h5(path):
    with h5py.File(path, "r") as h5:
        feats = h5["matrix/features"]
        ids = [x.decode() for x in feats["id"][:]]
        names = [x.decode() for x in feats["name"][:]]
        types = [x.decode() for x in feats["feature_type"][:]]
    keep = [i for i, t in enumerate(types) if t == "Gene Expression"]
    return {ids[i].split(".")[0] for i in keep}, {names[i] for i in keep}


def axis_pkc(path):
    with gzip.open(path, "rt") as fh:
        pkc = json.load(fh)
    sym = {t["DisplayName"] for t in pkc["Targets"] if t.get("CodeClass", "").startswith("Endogenous")}
    return set(), sym


def axis_csv_header(path, first_line_bytes=None):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    cols = [c.strip('"') for c in header]
    sym = {c for c in cols if c and c not in {"fov", "cell_ID", "cell_id"} and not c.startswith(CONTROL_PREFIXES)}
    controls = [c for c in cols if c.startswith(CONTROL_PREFIXES)]
    return set(), sym, controls


def axis_h5ad_var(path):
    with h5py.File(path, "r") as h5:
        var = h5["var"]
        key = var.attrs.get("_index", "_index")
        names = var[key][:]
    sym = {x.decode() if isinstance(x, bytes) else str(x) for x in names}
    return set(), {s for s in sym if not s.startswith(CONTROL_PREFIXES)}


def main(run_dir):
    run_dir = Path(run_dir)
    dl = run_dir / "downloads_untrusted"
    out = run_dir / "coverage"
    registry, members = load_programs()
    sym_to_ens, biotype = load_gencode()

    axes, notes = {}, {}
    axes["GSE298719_flex"] = axis_10x_features(dl / "gse298719/GSM9859944_features.tsv.gz")
    axes["GSE344087_citeseq"] = axis_10x_h5(dl / "gse344087/GSM9970284_S009.h5")
    axes["GSE255312_geomx_wta"] = axis_pkc(dl / "gse255312/GSE255312_Hs_R_NGS_WTA_v1.0.pkc.gz")
    for name, path in [("Fukumoto_cosmx_s1", dl / "fukumoto/expr_s1_head.bin"),
                       ("Fukumoto_cosmx_s2", dl / "fukumoto/expr_s2_head.bin"),
                       ("GSE292268_cosmx_slide1", dl / "gse292268/GSM8855042_Slide1_exprMat_file.csv.gz"),
                       ("GSE292268_cosmx_slide2", dl / "gse292268/GSM8855043_Slide2_exprMat_file.csv.gz")]:
        ens, sym, controls = axis_csv_header(path)
        axes[name] = (ens, sym)
        notes[name] = f"{len(sym)} targets; {len(controls)} control columns ({', '.join(sorted({c.rstrip('0123456789_') for c in controls})[:4])})"
    axes["Govaere_cosmx_reference"] = axis_h5ad_var(ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad")

    rows, missing_rows, summary = [], [], []
    for dataset, (ens_axis, sym_axis) in axes.items():
        n_obs = n_partial = n_untest = n_rep = 0
        for reg in registry:
            uid = reg["program_uid"]
            genes = members[uid]
            present, weight = [], 0.0
            for sym, w in genes:
                hit = sym in sym_axis or bool(sym_to_ens.get(sym, set()) & ens_axis)
                if hit:
                    present.append(sym)
                    weight += w
                elif uid in CONFIRMATORY:
                    missing_rows.append({"dataset": dataset, "program": CONFIRMATORY[uid], "gene": sym,
                                         "biotype_gencode_v49": biotype.get(sym, "NA"), "l1_weight": f"{w:.4f}"})
            if len(present) >= OBS_MIN_GENES and weight >= OBS_MIN_WEIGHT:
                state = "observable"; n_obs += 1
            elif present:
                state = "partial"; n_partial += 1
            else:
                state = "untestable"; n_untest += 1
            rep = len(present) >= REP_MIN_GENES and weight >= REP_MIN_WEIGHT
            n_rep += rep
            rows.append({"dataset": dataset, "program_uid": uid, "cell_type": reg["cell_type"],
                         "program_label": reg.get("module_name", ""), "n_confirmed_genes": len(genes),
                         "n_present": len(present), "presence_weight": f"{weight:.4f}",
                         "spatial_state_presence_only": state,
                         "replication_rule_presence_only": "pass" if rep else "fail"})
        summary.append({"dataset": dataset, "n_axis_genes": len(sym_axis | ens_axis) if not ens_axis else len(ens_axis),
                        "observable": n_obs, "partial": n_partial, "untestable": n_untest,
                        "pass_0.80_presence": n_rep, "note": notes.get(dataset, "")})

    def write(path, recs):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t")
            w.writeheader(); w.writerows(recs)

    write(out / "program_gene_coverage_all_datasets.tsv", rows)
    write(out / "confirmatory_missing_genes.tsv", missing_rows)
    write(out / "coverage_summary.tsv", summary)
    conf = [r for r in rows if r["program_uid"] in CONFIRMATORY]
    write(out / "confirmatory_program_coverage.tsv", conf)


if __name__ == "__main__":
    main(sys.argv[1])
