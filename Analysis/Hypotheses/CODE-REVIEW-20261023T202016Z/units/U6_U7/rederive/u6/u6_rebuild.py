#!/usr/bin/env python3
"""Code review U6: independent rebuild of the Gene Catalog r9d evidence classes from the adopted bulk fit
and COLOC r2, using the class rules in METHODS.md ("Joint testing universe and evidence classes") and
catalog_rules_r9.tsv. Not derived from build_catalog_r9*.py.

Rules used:
  RNA      supported  = in the 24,115-gene fit, BH q < 0.05 and |log2FC| > 0.5
           indeterminate = in the fit, not supported; untestable = no fit row
  genetics supported  = Tier-1/2 best SuSiE PP.H4 > 0.5
           untestable = no Tier-1/2 COLOC record, or susie_state_t12 == 'untestable'
           indeterminate = otherwise (ABF only or SuSiE tested at <= 0.5)
  class    both supported -> concordant; gen sup + RNA ind -> genetically_anchored;
           RNA sup + gen ind -> established_state_associated; both ind -> unresolved;
           either untestable -> untested
Identity: Ensembl version stripped; IDs listed in the r9d sidecar of IDs without a Catalog record are
removed (adjudication is an input, checked separately below); the six COLOC IDs that the Catalog carries
on another record are reported both ways.
Usage: u6_rebuild.py <out_dir>
"""
import sys, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
R9D = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out"
CAT = R9D / "catalog_r9d_gene_table.parquet"
SIDE = R9D / "sidecar_ids_without_catalog_record.tsv"
log = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
pd.DataFrame({"path": [str(p) for p in (BULK, COLOC, CAT, SIDE)], "sha256": [sha(p) for p in (BULK, COLOC, CAT, SIDE)]}).to_csv(OUT / "inputs_sha256.tsv", sep="\t", index=False)

# ---------------- bulk
b = pd.read_csv(BULK)
b["id"] = b.gene.str.split(".").str[0]
say("bulk rows", len(b), "unique versioned", b.gene.nunique(), "unique stripped", b.id.nunique())
assert b.id.is_unique, "version stripping merged IDs"
b["rna_sup"] = (b.padj < 0.05) & (b.logFC.abs() > 0.5)
say("bulk canonical calls", int(b.rna_sup.sum()), "NA padj", int(b.padj.isna().sum()))
# ---------------- COLOC
c = pd.read_csv(COLOC)
c["id"] = c.ensembl.astype(str).str.split(".").str[0]
say("COLOC tier12 rows", len(c), "unique stripped", c.id.nunique())
assert c.id.is_unique
c["gen_sup"] = c.coloc_best_susie_pp4 > 0.5
c["gen_unt"] = c.susie_state_t12 == "untestable"
say("COLOC SuSiE > 0.5:", int(c.gen_sup.sum()), "; >= 0.5:", int((c.coloc_best_susie_pp4 >= 0.5).sum()),
    "; state counts", c.susie_state_t12.value_counts().to_dict())
assert (c.gen_sup == (c.susie_state_t12 == "positive")).all()

# ---------------- Catalog and sidecar
W = pd.read_parquet(CAT)
W["id"] = W.ensembl_id.astype(str).str.split(".").str[0]
assert W.id.is_unique
S = pd.read_csv(SIDE, sep="\t", dtype=str, keep_default_na=False)
S["id"] = S.ensembl_id.str.split(".").str[0]
say("Catalog records", len(W), "; sidecar rows", len(S), S.source.value_counts().to_dict())

bulk_ids, coloc_ids, cat_ids = set(b.id), set(c.id), set(W.id)
side_bulk = set(S.id[S.source == "bulk_774_fit"]); side_coloc = set(S.id[S.source == "coloc_r2_tier12"])
side_prog = set(S.id[~S.source.isin(["bulk_774_fit", "coloc_r2_tier12"])])
U = bulk_ids | coloc_ids
say("union bulk|COLOC", len(U), "; bulk&COLOC (stable-ID join)", len(bulk_ids & coloc_ids))
say("sidecar bulk IDs in bulk fit", len(side_bulk & bulk_ids), "/", len(side_bulk), "; sidecar COLOC IDs in COLOC", len(side_coloc & coloc_ids), "/", len(side_coloc))
expected = U - side_bulk - side_coloc
say("union minus sidecar bulk/COLOC IDs", len(expected), "; Catalog", len(cat_ids))
only_exp = sorted(expected - cat_ids); only_cat = sorted(cat_ids - expected)
say("in union-minus-sidecar but not Catalog:", len(only_exp), only_exp[:20])
say("in Catalog but not union-minus-sidecar:", len(only_cat), only_cat[:20])
say("sidecar IDs that are also Catalog records:", len((side_bulk | side_coloc | side_prog) & cat_ids))

# six COLOC IDs carried on another record
oth = W[W.gen_estimate_from_other_ensembl_id.astype(str).str.startswith("ENSG")][["id", "symbol", "gen_estimate_from_other_ensembl_id"]]
say("records carrying genetics from another Ensembl ID:", len(oth)); say(oth.to_string(index=False))

# ---------------- independent states and classes
def classes(rekey):
    d = pd.DataFrame({"id": W.id})
    bi = b.set_index("id"); ci = c.set_index("id")
    gid = d.id.copy()
    if rekey:
        m = dict(zip(oth.id, oth.gen_estimate_from_other_ensembl_id.str.split(".").str[0]))
        gid = gid.map(lambda x: m.get(x, x))
    in_b = d.id.isin(bi.index); in_c = gid.isin(ci.index)
    rs = np.where(~in_b, "untestable", np.where(d.id.map(bi.rna_sup).fillna(False).astype(bool), "supported", "indeterminate"))
    gs_sup = gid.map(ci.gen_sup).fillna(False).astype(bool); gs_unt = gid.map(ci.gen_unt).fillna(False).astype(bool)
    gs = np.where(~in_c | gs_unt, "untestable", np.where(gs_sup, "supported", "indeterminate"))
    cl = np.select([(rs == "untestable") | (gs == "untestable"), (rs == "supported") & (gs == "supported"),
                    gs == "supported", rs == "supported"],
                   ["untested", "concordant", "genetically_anchored", "established_state_associated"], "unresolved")
    d["rna_state"] = rs; d["gen_state"] = gs; d["class"] = cl
    d["gen_unt_reason"] = np.where(~in_c, "no_tier12_coloc_record", np.where(gs_unt, "susie_untestable_in_ge1_tier12_study", ""))
    return d
for rekey in (True, False):
    d = classes(rekey)
    tag = "with six re-keyed COLOC IDs" if rekey else "strict same-stable-ID"
    say(f"--- {tag}")
    say("class counts", d["class"].value_counts().to_dict())
    tab = pd.crosstab(pd.Categorical(d.gen_state, ["supported", "indeterminate", "untestable"]),
                      pd.Categorical(d.rna_state, ["supported", "indeterminate", "untestable"]))
    say("3x3 gen rows x RNA cols\n" + tab.to_string())
    say("gen untestable reasons", d.gen_unt_reason.value_counts().to_dict())
    agree = (d["class"].values == W.primary_evidence_class.values)
    say("per-gene class agreement with r9d:", int(agree.sum()), "/", len(d))
    if not agree.all():
        dd = d[~agree].assign(r9d=W.primary_evidence_class.values[~agree], symbol=W.symbol.values[~agree])
        say(dd.to_string(index=False))
    # per-gene state agreement with r9d r8 state columns
    for mine, theirs in (("rna_state", "r8_rna_call_state"), ("gen_state", "r8_genetics_call_state")):
        th = W[theirs].astype(str).values
        say(f"{mine} vs {theirs}: agree", int((d[mine].values == th).sum()), "; r9d values", pd.Series(th).value_counts().to_dict())
    if rekey: D = d
D.to_csv(OUT / "u6_independent_classes.tsv.gz", sep="\t", index=False)

# ---------------- Fig 6D numbers
D["symbol"] = W.symbol.values
unt = D["class"] == "untested"
say("Fig6D untested split: no RNA call (not detected)", int((unt & (D.rna_state == "untestable")).sum()),
    "; no genetic call", int((unt & (D.gen_state == "untestable")).sum()))
say("hatched: untested with RNA supported", int((unt & (D.rna_state == "supported")).sum()), "; with genetics supported", int((unt & (D.gen_state == "supported")).sum()))
ci = c.set_index("id")
gidk = D.id.map(lambda x: dict(zip(oth.id, oth.gen_estimate_from_other_ensembl_id.str.split(".").str[0])).get(x, x))
abf = gidk.map(ci.coloc_best_abf_pp4); sus = gidk.map(ci.coloc_best_susie_pp4)
flag = (abf > 0.5) & sus.isna()
say("ABF > 0.5 and no Tier-1/2 SuSiE result, by class:", D[flag]["class"].value_counts().to_dict())
say("  same with SuSiE state == no_susie_result:", D[(abf > 0.5) & (gidk.map(ci.susie_state_t12) == "no_susie_result")]["class"].value_counts().to_dict())
say("1,854 calls split: in Catalog", int(D.rna_state.eq("supported").sum()), "; sidecar bulk IDs with a call",
    sorted(set(b.id[b.rna_sup]) & side_bulk))

# ---------------- sidecar groups: would any have changed a class?
bi = b.set_index("id")
rows = []
for _, r in S.iterrows():
    i = r.id
    rows.append({"source": r.source, "id": i, "reason": r.quarantine_reason, "symbol": r.symbol,
                 "in_bulk": i in bulk_ids, "rna_call": bool(bi.rna_sup.get(i, False)),
                 "in_coloc": i in coloc_ids, "susie_pp4": ci.coloc_best_susie_pp4.get(i, np.nan),
                 "susie_state": ci.susie_state_t12.get(i, ""), "abf_pp4": ci.coloc_best_abf_pp4.get(i, np.nan),
                 "symbol_has_catalog_record": r.symbol in set(W.symbol)})
SG = pd.DataFrame(rows)
SG.to_csv(OUT / "sidecar_groups_checked.tsv", sep="\t", index=False)
say("sidecar groups:\n" + SG.groupby(["source", "reason"]).agg(n=("id", "size"), in_bulk=("in_bulk", "sum"), rna_call=("rna_call", "sum"),
    in_coloc=("in_coloc", "sum"), susie_gt05=("susie_pp4", lambda x: int((x > 0.5).sum())), susie_tested=("susie_pp4", lambda x: int(x.notna().sum())),
    symbol_kept=("symbol_has_catalog_record", "sum")).to_string())
# second-ID rows: compare with the kept record of the same symbol
sec = SG[SG.reason.isin(["not_quarantined", ""]) | ~SG.source.isin(["bulk_774_fit", "coloc_r2_tier12"])]
say("non-quarantined sidecar rows by source:", sec.source.value_counts().to_dict())
say(SG[(SG.source != "program_membership") & SG.reason.eq("not_quarantined")].to_string(index=False))

# ---------------- universes
say("stable-ID join bulk & COLOC:", len(bulk_ids & coloc_ids))
say("Catalog joint (class != untested):", int((D["class"] != "untested").sum()))
j = (bulk_ids & coloc_ids)
say("stable-ID join minus SuSiE-untestable:", len(j - set(c.id[c.gen_unt])))
say("join IDs without Catalog record:", len(j - cat_ids), "; join IDs whose Catalog record is untested:",
    int(D[D.id.isin(j)]["class"].eq("untested").sum()))

# ---------------- missing coded as negative
W2 = W.copy()
say("r9d RNA call NA where no fit row:", int(W2.rna_canonical_call.isna().sum()), "; NA rows with RNA state untestable:",
    int((W2.rna_canonical_call.isna() & (D.rna_state == "untestable")).sum()))
for col in W2.columns:
    v = W2[col].astype(str).str.lower()
    hits = v.str.contains(r"tested_negative|^negative$|not_supported$", regex=True)
    if hits.any(): say("negative-looking values in", col, v[hits].value_counts().head(5).to_dict())
(OUT / "log.txt").write_text("\n".join(log) + "\n")
