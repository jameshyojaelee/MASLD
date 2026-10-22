#!/usr/bin/env python3
"""H71: Gene Catalog r8 evidence state of the Kumak 2026 and Elison 2025 gene lists (Part A),
and expression-stratified colocalization overlap of the Kumak lists in U_eval (Part B).
Prespec: PRESPEC_H71.md (sha256 b4509c7b...), copied into the run directory.

Order follows the prespec "Prespecified diagnostics": xlsx identity and hashes (S1), margins
(S3, G1), PC-A (S2), bins and G2 (S4), identity and role-swap checks, PC-B (S5), PC-C, NC-1
(S6), NC-2 (S2), permutation null means (S7 part), D3 formula check, then the primary T1
result, Holm over family F1, D1 to D9, and the reading.

Debug flags (never used for the reported run):
  --stop_before_outcome  exits after the diagnostics, before any O/E_strat, MH OR or T1 p.
  --debug_null           after PC-A, replaces each list's U_full members by random genes drawn
                         within the same AveExpr bins (U_eval: T1 bins, same L_k; rest of
                         U_full: U_full deciles). Outcome numbers are then not real.
"""
import argparse, hashlib, json, os, re, sys, time
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom, binomtest
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H71 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RR = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
R4 = RR + "round4-20261015T112148Z/"
PRESPEC_SHA = "b4509c7baba68fc09271b20a9934987ab51cf9713c00e36fef61f3c12bb673b6"
XLSX = "data/public_r4/biorxiv_2026.08.26.747338/media-1.xlsx"
KTSV = {"NAS": R4 + "gate/NG10/kumak_NAS.tsv", "Fibrosis": R4 + "gate/NG10/kumak_Fibrosis.tsv"}
ELISON = "data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv"
CAT = "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"
C = "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
A = "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
B = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
RDS = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds"
DISC = "docs/manuscript/draft/fig5_discussion.md"
COMP = "docs/competitor_analysis/kumak_konu_2026_masld_meta_network.md"
# (path, bytes, sha256): the 20 rows of the prespec input table
INPUTS = [
    (XLSX, 129503, "8f85caf1384db046e01a4d1ef0a34ac7c526a25c82a2af247fd1e53dae4b9c4e"),
    (KTSV["NAS"], 68063, "cb0798d717fa7fed155209915f73b6c3c4c857abc1b2d037f9527b7971c8cfbe"),
    (KTSV["Fibrosis"], 74167, "68f92d7a4e3ca9e3dc7777f0498601e0c657e85975191ea72ec081ac6e32460c"),
    (ELISON, 4668, "6457be4a0a800ad7b3eb3e771b89ff158164d6fa2b0ecfb3459593167c0dbf15"),
    (CAT + "passport_gene_index.parquet", 468456, "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b"),
    (CAT + "passport_evidence_long.parquet", 3194057, "a2d4131bd38464501e595a6554ae1b9956d9ae1fd503f57f240d87120d0c392b"),
    (C, 2245414, "6360984750864e7ab56593bce90e8ce665d781ec72a1ba52790d6ba53e83b8d6"),
    (A, 335138620, "acdfa306f1ee45383758bc2522602fc2828c4060080b7ad2d96134fc4a399887"),
    (B, 5049936, "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6"),
    (RDS, 33106829, "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b"),
    (R4 + "gate/NG10/coverage.py", 2252, "08d02480b501d1570002880070aed9a45a26478bec3cb2fdaff314cf7a3a335d"),
    (R4 + "gate/NG10/coverage.out", 2806, "a84087582edd705d660b85ab8648db092e91cdd4fdaee45c63c70506d49e33b4"),
    (DISC, 21030, "32e61db482c7df6b4e2504f20e7cdc2817ea4fc7cd0715c9e86802d562337777"),
    ("scripts/manuscript/build_external_scchromatin_annotation.R", 7961, "e3037b58123720d65159e4303d953aaba754173a6312992accbd9641c481594c"),
    (R4 + "write/H71/count_n.py", 3124, "5c211bd5e9c3cc4a2630f869e7ef8352d3ff038f5b6c817ab1d9bfe5a2d16d6a"),
    (R4 + "write/H71/count_n.out", 1703, "c6717a773bd65c00f18415d6d6f46c6bd28aa6e8a31bba18ee7f3c260e487dcf"),
    (R4 + "write/H71/mde.py", 1038, "7249ea67edbaca076196d3a76289bee75f34e61d5658c70502fa8329318b67a6"),
    (R4 + "write/H71/mde.out", 718, "a6084b38df3a80af9e6d7450452d06c3bdc77ceca98d786625ba94841feb069f"),
    (R4 + "write/H71/xlsx_check.py", 906, "60975080f9ec3fbef5c3a9de519f0fa67f7225dc91c7b9a8a12fb4664544bcc2"),
    (R4 + "write/H71/xlsx_check.out", 375, "5f42298757affdf34ddb595e9e1b47ba0a8fde574eed71d85139e927e7f06b40"),
]
LISTS = ["NAS", "Fibrosis"]
SHIFT = {"NAS": (15, -13), "Fibrosis": (14, -12)}          # PC-B planted shifts
G3 = {"NAS": dict(OR_up=1.742, OR_dn=0.485), "Fibrosis": dict(OR_up=1.747, OR_dn=0.481)}  # mde.out
NAMED = ["COL1A1", "COL3A1", "THBS2", "FBLN5", "PDGFA", "LAMA2", "LAMA3", "GNMT", "MAT1A", "CYP2C19"]
TOL = 1e-9
REL_DEP = "genetic columns from COLOC r2; rebuild on adopted COLOC rerun"

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20260923)
ap.add_argument("--nperm", type=int, default=10000)
ap.add_argument("--npcb", type=int, default=1000)
ap.add_argument("--nnc", type=int, default=1000)
ap.add_argument("--nnc2", type=int, default=1000)
ap.add_argument("--debug_null", action="store_true")
ap.add_argument("--stop_before_outcome", action="store_true")
args = ap.parse_args()
OUT = os.path.abspath(args.out)
os.chdir(ROOT)
DIAG = []
RES = {"debug_null": args.debug_null, "stop_before_outcome": args.stop_before_outcome,
       "nperm": args.nperm, "npcb": args.npcb, "nnc": args.nnc, "nnc2": args.nnc2, "seed": args.seed}
PARTA = {}            # filled after PC-A: table, summary
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def diag(section, key, value):
    DIAG.append((section, key, value))
    log(f"DIAG {section} | {key} = {value}")


def write_diag():
    pd.DataFrame(DIAG, columns=["section", "key", "value"]).to_csv(f"{OUT}/diagnostics.tsv", sep="\t", index=False)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def jdump(o):
    return json.dumps(o, default=lambda v: v.item() if hasattr(v, "item") else str(v))


def write_parta():
    if "table" in PARTA:
        PARTA["table"].to_csv(f"{OUT}/partA_requalification.tsv", sep="\t", index=False)
        PARTA["summary"].to_csv(f"{OUT}/partA_summary.tsv", sep="\t", index=False)


def finish(status, verdict, lines, write_a=True):
    RES["status"] = status
    RES["verdict"] = verdict
    write_diag()
    if write_a:
        write_parta()
    with open(f"{OUT}/results.json", "w") as f:
        json.dump(RES, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    hdr = ["# RESULT_H71", "", f"Run directory: `{OUT}`",
           f"Prespec: `{H71}/PRESPEC_H71.md` (sha256 {PRESPEC_SHA})",
           f"Debug null mode: {args.debug_null}. Stop before outcome: {args.stop_before_outcome}. "
           f"B: T1 permutations {args.nperm}, PC-B {args.npcb}, NC-1 {args.nnc}, NC-2 {args.nnc2}.", "",
           f"**Status:** {status}", "", f"**Verdict:** {verdict}", ""]
    with open(f"{OUT}/RESULT_H71.md", "w") as f:
        f.write("\n".join(hdr + lines) + "\n")
    log("STATUS", status, "|", verdict)
    sys.exit(0)


def stop(rule, msg, hold_parta=True):
    finish("stopped_by_rule", f"{rule}: {msg}. No number from this run is used"
           + ("; Part A is held." if hold_parta else "."),
           [f"Stop rule {rule} fired: {msg}", "", "See diagnostics.tsv for the values printed before the stop."],
           write_a=not hold_parta)


# ======================================================================= 1. xlsx identity, hashes (S1)
xc = open(f"{OUT}/xlsx_check.out").read().splitlines()
for s in LISTS:
    ln = [l for l in xc if l.startswith(f"{s} sheet")]
    ok = len(ln) == 1 and ln[0].rstrip().endswith("cell-identical True")
    diag("1_hash", f"xlsx_check {s}", ln[0] if ln else "MISSING")
    if not ok:
        stop("S1", f"xlsx check for {s} is not 'cell-identical True'")
for p in (f"{H71}/PRESPEC_H71.md", RR + "plans/PRESPEC_H71.md"):
    got = sha(p)
    diag("1_hash", p, got)
    if got != PRESPEC_SHA:
        stop("S1", f"prespec sha256 {got} != {PRESPEC_SHA} ({p})")
bad = []
for p, nb, h in INPUTS:
    g = sha(p); sz = os.path.getsize(p)
    diag("1_hash", p, "OK" if (g == h and sz == nb) else f"MISMATCH bytes {sz} sha256 {g}")
    if g != h or sz != nb:
        bad.append(p)
if bad:
    stop("S1", "sha256 or byte-count mismatch for " + "; ".join(bad))
RES["S1"] = "pass"

# ======================================================================= load
gi = pd.read_parquet(CAT + "passport_gene_index.parquet")
ev = pd.read_parquet(CAT + "passport_evidence_long.parquet")
tr = ev[ev.evidence_domain == "transcriptomics"].set_index("ensembl_id")
ge = ev[ev.evidence_domain == "genetics"].set_index("ensembl_id")
assert tr.index.is_unique and ge.index.is_unique
gi["g"] = gi.ensembl_id.astype(str).str.split(".").str[0]


def direction(k):
    return np.where(k.total_upregulated > k.total_downregulated, "up",
                    np.where(k.total_upregulated < k.total_downregulated, "down", "tie"))


kum = {}
for s in LISTS:
    k = pd.read_csv(KTSV[s], sep="\t")
    k["dir"] = direction(k)
    kum[s] = k
el = pd.read_csv(ELISON)

bulk = pd.read_csv(B)
bulk["g"] = bulk.gene.astype(str).str.split(".").str[0]
bulk = bulk.drop_duplicates("g")
bulk["deg"] = ((bulk.padj < 0.05) & (bulk.logFC.abs() > 0.5)).astype(bool)
coloc = pd.read_csv(C)
coloc["g"] = coloc.ensembl.astype(str).str.split(".").str[0]
diag("2_sizes", "coloc rows / unique g", f"{len(coloc)} / {coloc.g.nunique()}")
coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")

full = bulk.merge(coloc[["g", "coloc_best_susie_pp4", "susie_state_t12", "susie_driving_gwas"]], on="g", how="inner")
full["susie"] = pd.to_numeric(full.coloc_best_susie_pp4, errors="coerce")
full["gen"] = (full.susie >= 0.5).fillna(False).astype(bool)
full = full.sort_values("g").reset_index(drop=True)                    # sorted by Ensembl ID (prespec bins)
ueval = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)

# ======================================================================= 2. margins (S3) and G1
mism = []


def expect(key, got, want):
    ok = got == want
    diag("2_margins", key, f"{got}" + ("" if ok else f"  != count_n.out {want}"))
    if not ok:
        mism.append(key)


expect("catalog genes", len(gi), 28233)
expect("symbol_collision True", int(gi.symbol_collision.fillna(False).sum()), 0)
expect("dup symbols", int(gi.symbol.duplicated(keep=False).sum()), 0)
CN = {  # count_n.out
    "NAS": dict(rows=951, uniq=951, dir={"up": 654, "down": 297},
                GTEX={"Expressed in Liver": 853, "Highest in Liver": 91, "Not found": 6, "Not Expressed in Liver": 1},
                pub={True: 561, False: 390},
                unmatched=["C11orf80", "C12orf4", "CSF2RA", "GOLGA2P5", "LPAL2", "PLGLA", "TM4SF19-DYNLT2B"],
                matched=944, inbulk=938, infull=902, ineval=241,
                states={"no_susie_result": 661, "tested_le_0.5": 205, "positive": 36},
                q={0.0: -1.922, 0.1: 1.752, 0.5: 4.227, 0.9: 6.886, 1.0: 11.401}),
    "Fibrosis": dict(rows=1028, uniq=1028, dir={"up": 864, "down": 164},
                     GTEX={"Expressed in Liver": 924, "Highest in Liver": 84, "Not found": 18, "Not Expressed in Liver": 2},
                     pub={True: 607, False: 421},
                     unmatched=["ARPC4-TTLL3", "C11orf80", "C3P1", "COL4A2-AS1", "FKBP1A-SDCBP2", "FOLH1B", "FXYD6-FXYD2",
                                "IGLC1", "INMT-MINDY4", "KIAA1522", "MAGED4", "MAGED4B", "MIR198", "P2RX5-TAX1BP3", "ZNF204P"],
                     matched=1013, inbulk=1005, infull=943, ineval=228,
                     states={"no_susie_result": 714, "tested_le_0.5": 200, "positive": 28, "untestable": 1},
                     q={0.0: -1.952, 0.1: 0.604, 0.5: 2.978, 0.9: 6.978, 1.0: 12.256}),
}
LSET = {}            # matched Ensembl IDs per list (real)
MERGED = {}
for s in LISTS:
    k = kum[s]; t = CN[s]
    expect(f"{s} tsv rows", len(k), t["rows"]); expect(f"{s} unique genes", int(k.Gene.nunique()), t["uniq"])
    expect(f"{s} dir", k.dir.value_counts().to_dict(), t["dir"])
    expect(f"{s} GTEX", k.GTEX.value_counts().to_dict(), t["GTEX"])
    expect(f"{s} pubmed_mentioned", k.pubmed_mentioned.value_counts().to_dict(), t["pub"])
    m = k.merge(gi[["symbol", "g"]], left_on="Gene", right_on="symbol", how="left")
    MERGED[s] = m
    expect(f"{s} unmatched", sorted(m.loc[m.g.isna(), "Gene"].tolist()), t["unmatched"])
    LSET[s] = set(m.g.dropna())
expect("NAS & Fibrosis matched overlap", len(LSET["NAS"] & LSET["Fibrosis"]), 290)
expect("NAS | Fibrosis union", len(LSET["NAS"] | LSET["Fibrosis"]), 1667)
expect("U_full N", len(full), 14704); expect("U_full G", int(full.gen.sum()), 426)
expect("U_eval N", len(ueval), 3769); expect("U_eval G", int(ueval.gen.sum()), 426)
expect("U_eval genetic all positive", bool((ueval.gen == (ueval.susie_state_t12 == "positive")).all()), True)
expect("catalog genes in bulk table", int(gi.g.isin(bulk.g).sum()), 24001)
expect("catalog genes in U_full", int(gi.g.isin(full.g).sum()), 14683)
for s in LISTS:
    L = LSET[s]; t = CN[s]
    expect(f"{s} matched", len(L), t["matched"]); expect(f"{s} in bulk table", int(bulk.g.isin(L).sum()), t["inbulk"])
    expect(f"{s} in U_full", int(full.g.isin(L).sum()), t["infull"]); expect(f"{s} in U_eval", int(ueval.g.isin(L).sum()), t["ineval"])
    expect(f"{s} SuSiE state in U_full", full[full.g.isin(L)].susie_state_t12.value_counts().to_dict(), t["states"])
    expect(f"{s} AveExpr quantiles in U_eval", ueval[ueval.g.isin(L)].AveExpr.quantile([0, .1, .5, .9, 1]).round(3).to_dict(), t["q"])
expect("U_eval AveExpr quantiles", ueval.AveExpr.quantile([0, .1, .5, .9, 1]).round(3).to_dict(),
       {0.0: -3.684, 0.1: -0.03, 0.5: 3.728, 0.9: 6.372, 1.0: 13.575})
expect("Elison rows", len(el), 20)
expect("Elison types", el.elison_finding_type.value_counts().to_dict(),
       {"boundary_absent": 4, "cell_of_action": 4, "featured_target": 4, "grn_tf_masl": 3, "grn_tf_fibplus": 2,
        "grn_tf_decreased": 2, "colocalized_locus": 1})
if mism:
    stop("S3", "margins differ from count_n.out: " + "; ".join(mism))
g1 = all(int(ueval.g.isin(LSET[s]).sum()) >= 100 for s in LISTS) and len(ueval) == 3769 and int(ueval.gen.sum()) == 426
diag("5_gates", "G1 (each list >= 100 in U_eval; U_eval 3769 / 426)", "pass" if g1 else "FAIL")
RES["S3"] = "pass"; RES["G1"] = bool(g1)
if not g1:
    stop("S3", "G1 set sizes failed")

# ======================================================================= 3. PC-A reproduction (S2)
pca_fail = []


def pca(key, got, want):
    ok = got == want
    diag("3_PC-A", key, f"{got}" + ("" if ok else f"  != gate {want}"))
    if not ok:
        pca_fail.append(key)


# gate logic (coverage.py): symbol join after dropping duplicated symbols; Catalog transcriptomics estimate
g1x = gi.drop_duplicates("symbol", keep=False)
pca("catalog classes", gi.primary_evidence_class.value_counts().to_dict(),
    {"untested": 13551, "unresolved": 13127, "established_state_associated": 1129, "genetically_anchored": 400, "concordant": 26})
pca("symbols mapping >1 ensembl", int((gi.groupby("symbol").ensembl_id.nunique() > 1).sum()), 0)
trb = ev[ev.evidence_domain == "transcriptomics"][["ensembl_id", "estimate", "q_value", "call_state", "testability_state"]]
GATE = {
    "NAS": dict(matched=944, cls={"unresolved": 623, "established_state_associated": 243, "untested": 42,
                                  "genetically_anchored": 29, "concordant": 7, "NA": 7},
                n_dir=938, agree=931, n_q=927, agree_q=923, d1=258, ac=36,
                named={"COL1A1": ("up", "established_state_associated"), "COL3A1": ("up", "established_state_associated"),
                       "THBS2": ("up", "established_state_associated"), "FBLN5": ("up", "established_state_associated"),
                       "PDGFA": ("up", "established_state_associated"), "LAMA2": ("up", "established_state_associated"),
                       "LAMA3": ("up", "established_state_associated"), "GNMT": ("down", "established_state_associated"),
                       "MAT1A": ("down", "unresolved"), "CYP2C19": ("down", "established_state_associated")}),
    "Fibrosis": dict(matched=1013, cls={"unresolved": 651, "established_state_associated": 263, "untested": 71,
                                        "genetically_anchored": 20, "NA": 15, "concordant": 8},
                     n_dir=1005, agree=975, n_q=930, agree_q=922, d1=301, ac=28,
                     named={"COL1A1": ("up", "established_state_associated"), "COL3A1": ("up", "established_state_associated"),
                            "THBS2": ("up", "established_state_associated"), "FBLN5": ("up", "established_state_associated"),
                            "PDGFA": ("up", "established_state_associated"), "LAMA2": ("up", "established_state_associated"),
                            "MAT1A": ("down", "unresolved"), "CYP2C19": ("down", "established_state_associated")}),
}
SIGN = {}
for s in LISTS:
    k = kum[s]; t = GATE[s]
    m = k.merge(g1x, left_on="Gene", right_on="symbol", how="left")
    pca(f"{s} matched", int(m.ensembl_id.notna().sum()), t["matched"])
    pca(f"{s} classes", {("NA" if pd.isna(a) else a): int(b) for a, b in m.primary_evidence_class.value_counts(dropna=False).items()}, t["cls"])
    pca(f"{s} anchored+concordant", int(m.primary_evidence_class.isin(["genetically_anchored", "concordant"]).sum()), t["ac"])
    mb = m.merge(trb, on="ensembl_id", how="inner")
    tt = mb[(mb.dir != "tie") & mb.estimate.notna()].copy()
    agree = ((tt.dir == "up") & (tt.estimate > 0)) | ((tt.dir == "down") & (tt.estimate < 0))
    pca(f"{s} bulk tested with direction", len(tt), t["n_dir"]); pca(f"{s} sign agree", int(agree.sum()), t["agree"])
    bt = binomtest(int(agree.sum()), len(tt), 0.5, alternative="greater").pvalue
    diag("3_PC-A", f"{s} sign-agreement binomial p (nominal, one-sided)", f"{bt:.4g}")
    sig = tt[tt.q_value < 0.05]
    ag = ((sig.dir == "up") & (sig.estimate > 0)) | ((sig.dir == "down") & (sig.estimate < 0))
    pca(f"{s} q<0.05 subset", len(sig), t["n_q"]); pca(f"{s} q<0.05 agree", int(ag.sum()), t["agree_q"])
    for gene in NAMED:
        r = m[m.Gene == gene]
        want = t["named"].get(gene)
        got = (r.dir.iloc[0], r.primary_evidence_class.iloc[0]) if len(r) else None
        pca(f"{s} named {gene}", got, want)
    mm = m.dropna(subset=["ensembl_id"])
    pca(f"{s} D1 bulk supported among matched", int((tr.loc[mm.ensembl_id, "call_state"] == "supported").sum()), t["d1"])
    tt = tt.assign(agree=agree.values)
    SIGN[s] = tt[["ensembl_id", "dir", "estimate", "q_value", "agree"]].reset_index(drop=True)
pca("whole Catalog anchored+concordant", int(gi.primary_evidence_class.isin(["genetically_anchored", "concordant"]).sum()), 426)
me = el.merge(g1x, left_on="gene", right_on="symbol", how="left")
pca("Elison n / matched", (len(el), int(me.ensembl_id.notna().sum())), (20, 20))
pca("Elison classes", me.primary_evidence_class.value_counts().to_dict(),
    {"unresolved": 11, "genetically_anchored": 5, "established_state_associated": 2, "concordant": 2})
EL_GATE = {"RORA": "genetically_anchored", "HNF4A": "unresolved", "CUX2": "established_state_associated",
           "CREB5": "unresolved", "RELB": "established_state_associated", "THRB": "genetically_anchored",
           "ONECUT1": "unresolved", "EFHD1": "concordant", "KRT8": "concordant", "PPP1R3B": "genetically_anchored",
           "TRPS1": "genetically_anchored", "HLA-DQA1": "unresolved", "FCGR2B": "genetically_anchored",
           "FLACC1": "unresolved", "CEBPA": "unresolved", "XBP1": "unresolved", "PNPLA3": "unresolved",
           "TM6SF2": "unresolved", "HSD17B13": "unresolved", "MBOAT7": "unresolved"}
pca("Elison per-gene classes", dict(zip(me.gene, me.primary_evidence_class)), EL_GATE)
RES["PC_A"] = {"pass": not pca_fail, "failed": pca_fail}
if pca_fail:
    stop("S2", "PC-A reproduction failed: " + ", ".join(pca_fail))

# ======================================================================= Part A table (annotation, no inference)
coloc_i = coloc.set_index("g")
bulk_i = bulk.set_index("g")
full_g = set(full.g); ueval_g = set(ueval.g)
gi_sym = gi.set_index("symbol")


def gene_cols(sym):
    """Catalog, bulk and COLOC columns for one list symbol (exact symbol join, no alias rescue)."""
    r = {}
    if sym in gi_sym.index:
        row = gi_sym.loc[sym]; g = row["g"]
        r.update(ensembl_id=g, catalog_symbol=sym, match_state="matched", primary_evidence_class=row["primary_evidence_class"])
        if g in bulk_i.index:
            b = bulk_i.loc[g]
            r.update(bulk_log2fc=float(b.logFC), bulk_q=float(b.padj), bulk_canonical_deg=bool(b.deg), in_bulk_table=True)
        else:
            r.update(bulk_log2fc=np.nan, bulk_q=np.nan, bulk_canonical_deg=False, in_bulk_table=False)
        r.update(bulk_call_state=tr.loc[g, "call_state"], bulk_testability_state=tr.loc[g, "testability_state"],
                 bulk_testability_reason=tr.loc[g, "testability_reason"])
        if g in coloc_i.index:
            c = coloc_i.loc[g]
            r.update(coloc_best_susie_pp4=pd.to_numeric(c.coloc_best_susie_pp4, errors="coerce"), susie_state_t12=c.susie_state_t12)
        else:
            r.update(coloc_best_susie_pp4=np.nan, susie_state_t12="not_in_coloc_table")
        r.update(genetic_call_state=ge.loc[g, "call_state"], genetic_testability_state=ge.loc[g, "testability_state"],
                 genetic_testability_reason=ge.loc[g, "testability_reason"],
                 in_U_full=g in full_g, in_U_eval=g in ueval_g)
    else:
        r.update(ensembl_id=np.nan, catalog_symbol=np.nan, match_state="unmatched", primary_evidence_class=np.nan,
                 bulk_log2fc=np.nan, bulk_q=np.nan, bulk_canonical_deg=np.nan, in_bulk_table=False, bulk_call_state=np.nan,
                 bulk_testability_state=np.nan, bulk_testability_reason=np.nan, coloc_best_susie_pp4=np.nan,
                 susie_state_t12=np.nan, genetic_call_state=np.nan, genetic_testability_state=np.nan,
                 genetic_testability_reason=np.nan, in_U_full=False, in_U_eval=False)
    return r


rows = []
for s in LISTS:
    for x in kum[s].itertuples(index=False):
        r = dict(list=f"Kumak_{s}", list_symbol=x.Gene, kumak_direction=x.dir, elison_finding_type=np.nan,
                 elison_direction=np.nan, kumak_total=int(x.total), GTEX=x.GTEX, pubmed_mentioned=x.pubmed_mentioned)
        r.update(gene_cols(x.Gene))
        lf = r["bulk_log2fc"]
        r["sign_agrees_with_list"] = (np.nan if pd.isna(lf) else bool((x.dir == "up" and lf > 0) or (x.dir == "down" and lf < 0)))
        r["release_dependence"] = REL_DEP
        rows.append(r)
for x in el.itertuples(index=False):
    r = dict(list="Elison", list_symbol=x.gene, kumak_direction=np.nan, elison_finding_type=x.elison_finding_type,
             elison_direction=x.elison_direction, kumak_total=np.nan, GTEX=np.nan, pubmed_mentioned=np.nan)
    r.update(gene_cols(x.gene))
    r["sign_agrees_with_list"] = np.nan          # Elison directions are GRN activity or chromatin, not RNA level
    r["release_dependence"] = REL_DEP
    rows.append(r)
TAB = pd.DataFrame(rows)
COLS = ["list", "list_symbol", "kumak_direction", "elison_finding_type", "elison_direction", "kumak_total", "GTEX",
        "pubmed_mentioned", "ensembl_id", "catalog_symbol", "match_state", "primary_evidence_class", "bulk_log2fc", "bulk_q",
        "bulk_call_state", "bulk_testability_state", "bulk_testability_reason", "bulk_canonical_deg", "in_bulk_table",
        "sign_agrees_with_list", "coloc_best_susie_pp4", "susie_state_t12", "genetic_call_state", "genetic_testability_state",
        "genetic_testability_reason", "in_U_full", "in_U_eval", "release_dependence"]
TAB = TAB[COLS]
# Part A consistency with the gate join
for s in LISTS:
    t = TAB[TAB.list == f"Kumak_{s}"]
    diag("3_PartA", f"{s} rows / matched / in U_eval", f"{len(t)} / {int((t.match_state == 'matched').sum())} / {int(t.in_U_eval.sum())}")
    diag("3_PartA", f"{s} max |bulk_log2fc - Catalog estimate|",
         f"{float((t.dropna(subset=['ensembl_id']).set_index('ensembl_id').bulk_log2fc - tr.estimate.astype(float)).abs().max()):.3g}")

# summary rows and D1 (never tested)
SUM = []
D1 = {}
for lst in ["Kumak_NAS", "Kumak_Fibrosis", "Elison"]:
    t = TAB[(TAB.list == lst) & (TAB.match_state == "matched")]
    n = len(t)
    sup = t.bulk_call_state == "supported"
    SUM.append((lst, "list genes", "", int((TAB.list == lst).sum())))
    SUM.append((lst, "matched", "", n))
    SUM.append((lst, "bulk supported (D1)", "count", int(sup.sum())))
    SUM.append((lst, "bulk supported (D1)", "fraction of matched", round(float(sup.mean()), 4)))
    for c, v in t.primary_evidence_class.value_counts().items():
        SUM.append((lst, "primary class", c, int(v)))
    for c, v in t.susie_state_t12.value_counts().items():
        SUM.append((lst, "susie_state_t12 (all matched)", c, int(v)))
    for c, v in t[t.in_U_full].susie_state_t12.value_counts().items():
        SUM.append((lst, "susie_state_t12 (in U_full)", c, int(v)))
    ac = t[t.primary_evidence_class.isin(["genetically_anchored", "concordant"])]
    for c, v in ac.susie_state_t12.value_counts().items():
        SUM.append((lst, "anchored+concordant by susie_state_t12", c, int(v)))
    if lst != "Elison":
        q05_small = (t.bulk_q < 0.05) & (t.bulk_log2fc.abs() <= 0.5)
        D1[lst] = dict(matched=n, supported=int(sup.sum()), frac=float(sup.mean()),
                       by_class=t[sup].primary_evidence_class.value_counts().to_dict(),
                       by_direction=t[sup].kumak_direction.value_counts().to_dict(),
                       direction_totals=t.kumak_direction.value_counts().to_dict(),
                       q05_small_effect=int(q05_small.sum()),
                       q05_small_effect_by_class=t[q05_small].primary_evidence_class.value_counts().to_dict())
        for c, v in D1[lst]["by_class"].items():
            SUM.append((lst, "D1 supported by class", c, int(v)))
        for c, v in D1[lst]["by_direction"].items():
            SUM.append((lst, "D1 supported by Kumak direction", c, int(v)))
        SUM.append((lst, "D1 q<0.05 and |log2FC|<=0.5", "count", D1[lst]["q05_small_effect"]))
        diag("3_PartA", f"D1 {lst}", jdump(D1[lst]))
PARTA["table"] = TAB
PARTA["summary"] = pd.DataFrame(SUM, columns=["list", "quantity", "level", "value"])
RES["D1"] = D1


# ======================================================================= helpers
def make_bins(df, col, k):
    """df must already be sorted by Ensembl ID; ties in col are broken by that order."""
    return pd.qcut(df[col].rank(method="first"), k, labels=False).to_numpy()


def bin_margins(bins, gen, lst):
    ks = np.unique(bins)
    n = np.array([(bins == k).sum() for k in ks])
    Gk = np.array([gen[bins == k].sum() for k in ks])
    Lk = np.array([lst[bins == k].sum() for k in ks])
    Ok = np.array([(gen & lst)[bins == k].sum() for k in ks])
    return n, Gk, Lk, Ok


def exact_null(n, Gk, Lk):
    """pmf of the total overlap: convolution of per-bin hypergeometrics (population n_k, L_k successes, G_k draws)."""
    pmf = np.array([1.0])
    for nk, gk, lk in zip(n, Gk, Lk):
        x = np.arange(0, min(gk, lk) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, nk, lk, gk))
    return pmf


def exact_p(pmf, E, O):
    x = np.arange(len(pmf))
    return float(min(1.0, pmf[np.abs(x - E) >= abs(O - E) - TOL].sum()))


def mh(n, Gk, Lk, Ok):
    tabs = [np.array([[o, g - o], [l - o, nn - g - l + o]], dtype=float) for nn, g, l, o in zip(n, Gk, Lk, Ok)]
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint(alpha=0.05)
    r = dict(MH_OR=float(st.oddsratio_pooled), MH_lo=float(lo), MH_hi=float(hi))
    try:
        bd = st.test_equal_odds(adjust=False); r.update(BD_chi2=float(bd.statistic), BD_p=float(bd.pvalue))
    except Exception as e:  # noqa: BLE001
        r.update(BD_chi2=np.nan, BD_p=np.nan, BD_note=str(e))
    try:
        cmh = st.test_null_odds(correction=False); r.update(CMH_chi2=float(cmh.statistic), CMH_p=float(cmh.pvalue))
    except Exception as e:  # noqa: BLE001
        r.update(CMH_chi2=np.nan, CMH_p=np.nan)
    return r


def strat_test(bins, gen, lst, with_mh=True):
    n, Gk, Lk, Ok = bin_margins(bins, gen, lst)
    keep = n > 0
    n, Gk, Lk, Ok = n[keep], Gk[keep], Lk[keep], Ok[keep]
    E = float((Gk * Lk / n).sum())
    O = int(Ok.sum())
    pmf = exact_null(n, Gk, Lk)
    r = dict(n_strata=len(n), N=int(n.sum()), G=int(Gk.sum()), L=int(Lk.sum()), O=O, E_strat=E,
             O_over_E=O / E if E > 0 else np.nan, p_exact=exact_p(pmf, E, O),
             null_mean=float((np.arange(len(pmf)) * pmf).sum()))
    if with_mh:
        r.update(mh(n, Gk, Lk, Ok))
    return r


def perm_null(bins, gen, lst, nperm, rng):
    idx = [np.where(bins == k)[0] for k in np.unique(bins)]
    sims = np.empty(nperm, dtype=int)
    for r in range(nperm):
        tot = 0
        for ii in idx:
            tot += int((rng.permutation(gen[ii]) & lst[ii]).sum())
        sims[r] = tot
    return sims


def perm_p(sims, E, O):
    return (int((np.abs(sims - E) >= abs(O - E) - TOL).sum()) + 1) / (len(sims) + 1)


def unstrat(gen, lst):
    N = len(gen); G = int(gen.sum()); Lt = int(lst.sum()); O = int((gen & lst).sum())
    orr, p = stats.fisher_exact([[O, G - O], [Lt - O, N - G - Lt + O]])
    return dict(N=N, G=G, L=Lt, O=O, E0=G * Lt / N, O_over_E0=O / (G * Lt / N), fisher_OR=float(orr), fisher_p=float(p))


ss = np.random.SeedSequence(args.seed).spawn(6)
RNG = {nm: np.random.default_rng(s) for nm, s in zip(
    ["T1_NAS", "T1_Fibrosis", "PCB", "NC1_NAS", "NC1_Fibrosis", "NC2"], ss)}

# ======================================================================= debug null lists (never in the reported run)
if args.debug_null:
    rng_dbg = np.random.default_rng(999)
    bins_dbg_e = make_bins(ueval, "AveExpr", 10)
    rest = full[~full.g.isin(ueval_g)].reset_index(drop=True)
    bins_dbg_r = make_bins(rest, "AveExpr", 10)
    for s in LISTS:
        L = LSET[s]
        inE = ueval.g.isin(L).to_numpy(); inR = rest.g.isin(L).to_numpy()
        new = set(g for g in L if g not in full_g)
        for k in np.unique(bins_dbg_e):
            ii = np.where(bins_dbg_e == k)[0]
            new |= set(ueval.g.to_numpy()[rng_dbg.choice(ii, int(inE[ii].sum()), replace=False)])
        for k in np.unique(bins_dbg_r):
            ii = np.where(bins_dbg_r == k)[0]
            new |= set(rest.g.to_numpy()[rng_dbg.choice(ii, int(inR[ii].sum()), replace=False)])
        LSET[s] = new
    log("DEBUG NULL: list memberships in U_full replaced by within-bin random genes. Part B numbers below are not real.")

gen_e = ueval.gen.to_numpy(); gen_f = full.gen.to_numpy(); deg_f = full.deg.to_numpy(); deg_e = ueval.deg.to_numpy()
inE = {s: ueval.g.isin(LSET[s]).to_numpy() for s in LISTS}
inF = {s: full.g.isin(LSET[s]).to_numpy() for s in LISTS}

# ======================================================================= 4. bins and G2 (S4)
bins10 = make_bins(ueval, "AveExpr", 10)
bt = []
for k in np.unique(bins10):
    ii = bins10 == k
    bt.append(dict(bin=int(k), AveExpr_min=float(ueval.AveExpr[ii].min()), AveExpr_max=float(ueval.AveExpr[ii].max()),
                   n_k=int(ii.sum()), G_k=int(gen_e[ii].sum()), L_k_NAS=int(inE["NAS"][ii].sum()),
                   L_k_Fibrosis=int(inE["Fibrosis"][ii].sum())))
bt = pd.DataFrame(bt)
bt.to_csv(f"{OUT}/bins_aveexpr10_Ueval.tsv", sep="\t", index=False)
for r in bt.itertuples():
    diag("4_bins", f"U_eval bin {r.bin}", f"AveExpr [{r.AveExpr_min:.4f}, {r.AveExpr_max:.4f}] n {r.n_k} G {r.G_k} "
         f"L_NAS {r.L_k_NAS} L_Fib {r.L_k_Fibrosis}")
diag("4_bins", "AveExpr ties in U_eval", int(ueval.AveExpr.duplicated().sum()))
LS = {}
for s in LISTS:
    _, Gk, Lk, _ = bin_margins(bins10, gen_e, inE[s])
    ok10 = int(((Gk >= 5) & (Lk >= 5)).sum())
    K = None; ok5 = None
    if ok10 >= 7:
        K = 10
    else:
        b5 = make_bins(ueval, "AveExpr", 5)
        _, G5, L5, _ = bin_margins(b5, gen_e, inE[s])
        ok5 = int(((G5 >= 5) & (L5 >= 5)).sum())
        K = 5 if ok5 >= 4 else None
    diag("5_gates", f"G2 {s}: 10-bin bins meeting rule / 5-bin / K", f"{ok10} / {ok5} / {K}")
    LS[s] = dict(K=K, ok10=ok10, ok5=ok5, untestable=None)
    if K is None:
        LS[s]["untestable"] = "S4"
        continue
    bins = bins10 if K == 10 else make_bins(ueval, "AveExpr", K)
    n, Gk, Lk, _ = bin_margins(bins, gen_e, inE[s])
    LS[s].update(bins=bins, n=n, Gk=Gk, Lk=Lk, pmf=exact_null(n, Gk, Lk), E=float((Gk * Lk / n).sum()))
RES["G2"] = {s: dict(K=LS[s]["K"], ok10=LS[s]["ok10"], ok5=LS[s]["ok5"]) for s in LISTS}

# ======================================================================= 5/6. identity and role swap
for s in LISTS:
    if LS[s]["untestable"]:
        continue
    one = np.zeros(len(ueval), dtype=int)
    n1, G1_, L1_, _ = bin_margins(one, gen_e, inE[s])
    idd = float((G1_ * L1_ / n1).sum()) - G1_.sum() * L1_.sum() / len(ueval)
    pmf = LS[s]["pmf"]; pmf_sw = exact_null(LS[s]["n"], LS[s]["Lk"], LS[s]["Gk"])
    m1 = float((np.arange(len(pmf)) * pmf).sum()); m2 = float((np.arange(len(pmf_sw)) * pmf_sw).sum())
    Lmax = max(len(pmf), len(pmf_sw)); pad = lambda v: np.pad(v, (0, Lmax - len(v)))
    diag("6_identity", f"{s} one-bin E_strat - G*L/N", f"{idd:.3g}")
    diag("7_roleswap", f"{s} null mean diff / max pmf diff", f"{m1 - m2:.3g} / {np.abs(pad(pmf) - pad(pmf_sw)).max():.3g}")
    diag("7_roleswap", f"{s} exact null mean - analytic E_strat", f"{m1 - LS[s]['E']:.3g}")
    LS[s].update(identity_ok=abs(idd) < 1e-9, roleswap_ok=abs(m1 - m2) < 1e-9)

# ======================================================================= 7. PC-B (S5): planted shifts, stream 2
pcb_rows = []
for s in LISTS:
    for sgn in (0, 1):
        sh = SHIFT[s][sgn]
        lab = f"{s} {'up' if sh > 0 else 'down'} {sh:+d}"
        if LS[s]["untestable"]:
            continue
        bins = LS[s]["bins"]; lst = inE[s]; pmf = LS[s]["pmf"]; E = LS[s]["E"]
        P_LOOK = np.array([exact_p(pmf, E, xx) for xx in range(len(pmf))])
        idx = [np.where(bins == k)[0] for k in np.unique(bins)]
        det = 0; det_side = 0; failed = 0
        for rep in range(args.npcb):
            gen = np.zeros(len(gen_e), dtype=bool)
            for ii in idx:
                gen[ii] = RNG["PCB"].permutation(gen_e[ii])
            ok = True
            for _ in range(abs(sh)):
                if sh > 0:   # genetic non-list <-> non-genetic list
                    A_ = [ii[gen[ii] & ~lst[ii]] for ii in idx]; B_ = [ii[~gen[ii] & lst[ii]] for ii in idx]
                else:        # genetic list <-> non-genetic non-list
                    A_ = [ii[gen[ii] & lst[ii]] for ii in idx]; B_ = [ii[~gen[ii] & ~lst[ii]] for ii in idx]
                w = np.array([len(a) * len(b) for a, b in zip(A_, B_)], dtype=float)
                if w.sum() == 0:
                    ok = False; break
                kk = RNG["PCB"].choice(len(idx), p=w / w.sum())
                a = RNG["PCB"].choice(A_[kk]); b = RNG["PCB"].choice(B_[kk])
                gen[a] = False; gen[b] = True
            if not ok:
                failed += 1; pcb_rows.append((lab, rep, -1, np.nan, False)); continue
            tot = int((gen & lst).sum()); p = P_LOOK[tot]; d = bool(p < 0.025)
            det += d; det_side += int(d and np.sign(tot - E) == np.sign(sh))
            pcb_rows.append((lab, rep, tot, p, d))
        rate = det / args.npcb
        LS[s][f"pcb_{'up' if sh > 0 else 'down'}"] = dict(shift=sh, rate=rate, rate_planted_side=det_side / args.npcb, failed=failed)
        band = "met" if 0.60 <= rate <= 0.85 else ("above (reported, not a failure)" if rate > 0.85 else
                                                   ("missed prediction (0.50 to 0.60)" if rate >= 0.50 else "below 0.50 (S5)"))
        diag("8_PC-B", lab, f"detection rate {rate:.3f} (planted side {det_side / args.npcb:.3f}), failed swaps {failed}; prediction 0.60-0.85: {band}")
pd.DataFrame(pcb_rows, columns=["case", "rep", "total", "p_exact", "detected"]).to_csv(f"{OUT}/pcb_replicates.tsv", sep="\t", index=False)
for s in LISTS:
    if LS[s]["untestable"]:
        continue
    ok = LS[s]["pcb_up"]["rate"] >= 0.5 and LS[s]["pcb_down"]["rate"] >= 0.5
    LS[s]["pcb_pass"] = bool(ok)
    diag("8_PC-B", f"{s} S5 triggered", not ok)
    if not ok:
        LS[s]["untestable"] = "S5"

# ======================================================================= 7. PC-C: list vs canonical DEG in U_full
binsF = make_bins(full, "AveExpr", 10)
for s in LISTS:
    pcc = strat_test(binsF, deg_f, inF[s])
    ok = bool(pcc["O_over_E"] > 2.0 and pcc["p_exact"] < 1e-10)
    LS[s]["pcc"] = pcc; LS[s]["pcc_pass"] = ok
    diag("8_PC-C", f"{s} U_full N / DEG / list / O / E_strat / O/E / exact p / MH OR",
         f"{pcc['N']} / {pcc['G']} / {pcc['L']} / {pcc['O']} / {pcc['E_strat']:.3f} / {pcc['O_over_E']:.3f} / {pcc['p_exact']:.3g} / {pcc['MH_OR']:.3f}")
    diag("8_PC-C", f"{s} pass (O/E > 2.0 and p < 1e-10)", ok)


# ======================================================================= 7. NC-1: expression-matched random lists (S6)
def nc1(s, rng):
    bins = LS[s]["bins"]; lst = inE[s]
    ae = ueval.AveExpr.to_numpy(); gids = ueval.g.to_numpy()
    li = np.where(lst)[0]
    pool = np.where(~lst)[0]
    pool = pool[np.argsort(ae[pool], kind="stable")]
    pv = ae[pool]; M = len(pool)
    keep = ~lst                                      # scoring universe: U_eval minus the Kumak list
    b2 = bins[keep]; g2 = gen_e[keep]
    G2 = int(g2.sum()); N2 = int(keep.sum())
    pos_in_keep = np.cumsum(keep) - 1                # index of a U_eval gene inside the reduced universe
    rows = []
    for _ in range(args.nnc):
        used = np.zeros(M, dtype=bool)
        chosen = []; gaps = []
        for gidx in rng.permutation(li):
            v = ae[gidx]
            pos = int(np.searchsorted(pv, v))
            r = pos
            while r < M and used[r]:
                r += 1
            l = pos - 1
            while l >= 0 and used[l]:
                l -= 1
            dr = pv[r] - v if r < M else np.inf
            dl = v - pv[l] if l >= 0 else np.inf
            if dl < dr:
                c = l
            elif dr < dl:
                c = r
            else:                                     # equal distance: lower Ensembl ID
                c = l if gids[pool[l]] < gids[pool[r]] else r
            used[c] = True; chosen.append(pool[c]); gaps.append(min(dl, dr))
        chosen = np.array(chosen)
        R = np.zeros(N2, dtype=bool); R[pos_in_keep[chosen]] = True
        a1 = ae[chosen]; a0 = ae[li]
        smd = (a1.mean() - a0.mean()) / np.sqrt((a1.var(ddof=1) + a0.var(ddof=1)) / 2)
        t = strat_test(b2, g2, R, with_mh=False)
        rows.append(dict(max_gap=float(max(gaps)), smd=float(smd), O=t["O"], E_strat=t["E_strat"], strat_OE=t["O_over_E"],
                         unstrat_OE=t["O"] / (G2 * len(chosen) / N2), p_exact=t["p_exact"]))
    return pd.DataFrame(rows)


nc_all = []
for s in LISTS:
    if LS[s]["untestable"]:
        continue
    t = nc1(s, RNG[f"NC1_{s}"])
    t.insert(0, "list", s); nc_all.append(t)
    r = dict(mean_strat_OE=float(t.strat_OE.mean()), fpr=float((t.p_exact < 0.05).mean()),
             mean_unstrat_OE=float(t.unstrat_OE.mean()), max_gap_median=float(t.max_gap.median()),
             max_gap_max=float(t.max_gap.max()), smd_mean=float(t.smd.mean()), smd_maxabs=float(t.smd.abs().max()),
             n_distinct_O=int(t.O.nunique()), O_min=int(t.O.min()), O_max=int(t.O.max()))
    pred = bool(0.95 <= r["mean_strat_OE"] <= 1.05 and 0.02 <= r["fpr"] <= 0.08)
    s6 = bool(r["fpr"] > 0.10 or not (0.90 <= r["mean_strat_OE"] <= 1.10))
    LS[s].update(nc1=r, nc1_pred_met=pred, S6=s6)
    diag("8_NC-1", s, jdump({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
    diag("8_NC-1", f"{s} prediction met (mean O/E 0.95-1.05, fpr 0.02-0.08)", pred)
    diag("8_NC-1", f"{s} S6 triggered (fpr > 0.10 or mean O/E outside 0.90-1.10)", s6)
if nc_all:
    pd.concat(nc_all).to_csv(f"{OUT}/nc1_lists.tsv", sep="\t", index=False)

# ======================================================================= 7. NC-2: permuted direction labels (S2), stream 5
nc2_rows = []
nc2_fail = []
for s in LISTS:
    t = SIGN[s]
    up = (t.dir == "up").to_numpy(); pos_ = (t.estimate > 0).to_numpy()
    obs = float((up == pos_).mean())
    p_up = up.mean(); f_pos = pos_.mean()
    exp_ = p_up * f_pos + (1 - p_up) * (1 - f_pos)
    perms = np.array([float((RNG["NC2"].permutation(up) == pos_).mean()) for _ in range(args.nnc2)])
    nc2_rows += [(s, i, v) for i, v in enumerate(perms)]
    pred = bool(abs(perms.mean() - exp_) <= 0.02)
    above = bool(obs > perms.max())
    kfull = kum[s]
    p_up_list = float((kfull.dir == "up").mean())
    diag("8_NC-2", s, f"n {len(t)} observed {obs:.4f} ({int((up == pos_).sum())}/{len(t)}); perm mean {perms.mean():.4f} "
         f"max {perms.max():.4f}; expected p_up*f_pos+p_down*f_neg {exp_:.4f} (p_up {p_up:.4f} among bulk-matched, "
         f"{p_up_list:.4f} over the whole list; f_pos {f_pos:.4f}); prediction within 0.02: {pred}; observed above every permuted value: {above}")
    LS[s]["nc2"] = dict(observed=obs, perm_mean=float(perms.mean()), perm_max=float(perms.max()), expected=exp_,
                        pred_met=pred, observed_above_all=above)
    if not above:
        nc2_fail.append(s)
pd.DataFrame(nc2_rows, columns=["list", "perm", "agreement"]).to_csv(f"{OUT}/nc2_perms.tsv", sep="\t", index=False)
if nc2_fail:
    stop("S2", "NC-2: observed sign agreement not above every permuted value for " + ", ".join(nc2_fail))

# ======================================================================= 8. T1 permutation null means (S7 part), streams 0 and 1
SIMS = {}
for s in LISTS:
    if LS[s]["untestable"]:
        continue
    SIMS[s] = perm_null(LS[s]["bins"], gen_e, inE[s], args.nperm, RNG[f"T1_{s}"])
    diag("9_T1perm", f"{s} permutation mean - analytic E_strat", f"{SIMS[s].mean() - LS[s]['E']:.4f}")
if SIMS:
    pd.DataFrame(SIMS).to_csv(f"{OUT}/t1_perm_sims.tsv", sep="\t", index=False)

# ======================================================================= 9. D3 formula check; D6 cluster count (margins)
ce = pd.read_csv(f"{OUT}/control_expression.tsv", sep="\t")
chk = bulk[["gene", "AveExpr"]].merge(ce, on="gene", how="inner")
d3_diff = float(np.abs(chk.ave_all_units - chk.AveExpr).max())
d3_ok = bool(len(chk) == 24115 and d3_diff < 1e-8)
diag("10_D3", "genes matched / max |all-unit mean - AveExpr|", f"{len(chk)} / {d3_diff:.3g}")
diag("10_D3", "formula check pass (< 1e-8)", d3_ok)
diag("10_D3", "Spearman control-only vs AveExpr (24,115 genes)", f"{stats.spearmanr(chk.ave_control, chk.AveExpr)[0]:.4f}")
RES["D3_formula_ok"] = d3_ok

allg = pd.read_csv(A, usecols=["gwas_name", "ensembl", "chr", "PP.H4.susie", "top_snp"], dtype={"chr": str, "top_snp": str})
allg["g"] = allg.ensembl.astype(str).str.split(".").str[0]
gg = ueval[ueval.gen][["g", "susie", "susie_driving_gwas"]]
pos = gg.merge(allg, left_on=["g", "susie_driving_gwas"], right_on=["g", "gwas_name"], how="left")
del allg
dup = pos.g.duplicated(keep=False)
n_dup = int(pos.g[dup].nunique())
if n_dup:
    pos = pos[~dup | ((pos["PP.H4.susie"] - pos.susie).abs() < 1e-9)].sort_values(["g", "top_snp"]).drop_duplicates("g")
diag("10_D6", "genetic genes / joined rows with top_snp / genes with duplicate rows",
     f"{len(gg)} / {int(pos.top_snp.notna().sum())} / {n_dup}")
diag("10_D6", "max |PP.H4.susie - coloc_best_susie_pp4|", f"{float((pos['PP.H4.susie'] - pos.susie).abs().max()):.3g}")
sp = pos.top_snp.str.split(":", expand=True)
pos["c"] = sp[0]; pos["p"] = pd.to_numeric(sp[1])
diag("10_D6", "top_snp chromosome equals chr column", f"{int((pos.c == pos['chr'].astype(str)).sum())} of {len(pos)}")
lab = pd.Series(index=pos.g.values, dtype=object)
for ch, grp in pos.sort_values(["c", "p", "g"]).groupby("c"):
    cid = np.concatenate([[0], np.cumsum(np.diff(grp.p.to_numpy()) > 500_000)])
    lab[grp.g.values] = [f"{ch}_{i}" for i in cid]
pos["cluster"] = lab[pos.g].values
nos = pos.cluster.isna()
pos.loc[nos, "cluster"] = "nosnp_" + pos.loc[nos, "g"]      # a gene without top_snp is its own cluster
diag("10_D6", "genetic genes without top_snp (own cluster)", int(nos.sum()))
diag("10_D6", "clusters among the 426 genetic genes (single linkage, 500 kb)", int(pos.cluster.nunique()))

write_diag()
log("Diagnostics written; entering outcome section.")
if args.stop_before_outcome:
    finish("debug_stop_before_outcome", "Debug run stopped before the outcome section; no outcome computed.", [])

# ======================================================================= OUTCOME: T1 per list
T1 = {}
for s in LISTS:
    if LS[s]["untestable"]:
        continue
    r = strat_test(LS[s]["bins"], gen_e, inE[s])
    sims = SIMS[s]
    pp = perm_p(sims, r["E_strat"], r["O"])
    band = 3 * np.sqrt(r["p_exact"] * (1 - r["p_exact"]) / args.nperm) + 1 / (args.nperm + 1)
    md = float(sims.mean() - r["E_strat"])
    r.update(K=LS[s]["K"], p_perm=pp, perm_band=band, perm_mean_minus_E=md,
             S7=bool(abs(pp - r["p_exact"]) > band or abs(md) >= 0.1),
             rbg_excludes_1=bool(r["MH_lo"] > 1 or r["MH_hi"] < 1))
    T1[s] = r
    log(f"T1 {s}", jdump({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}))
RES["T1"] = T1
s7 = [s for s in T1 if T1[s]["S7"]]
if s7:
    RES["S7"] = s7
    lines = ["Stop rule S7 fired for " + ", ".join(s7) + ". No Part B number is used (T1, D2 to D6 and the Part B controls).",
             "Part A and D1 are reported below.", ""]
    for lst, d in D1.items():
        lines.append(f"- D1 {lst}: {d['supported']} of {d['matched']} matched genes bulk-supported ({100 * d['frac']:.1f}%).")
    finish("stopped_by_rule", "S7: permutation check of T1 disagrees with the exact p. Part A and D1 are still reported.", lines)

# Holm over family F1 (two tests); an untestable list enters as p = 1 (fixed before the outcome; DEVIATIONS.md)
p_raw = [T1[s]["p_exact"] if s in T1 else 1.0 for s in LISTS]
_, p_holm, _, _ = multipletests(p_raw, alpha=0.05, method="holm")
for s, ph in zip(LISTS, p_holm):
    if s in T1:
        T1[s]["p_holm"] = float(ph)

# ======================================================================= descriptives D2 to D9 (nominal)
D = {}
D["D2"] = {s: unstrat(gen_e, inE[s]) for s in LISTS}
D["D2_catalog_wide"] = {s: dict(anchored_concordant=GATE[s]["ac"], matched=GATE[s]["matched"],
                                frac=GATE[s]["ac"] / GATE[s]["matched"]) for s in LISTS}
D["D2_catalog_wide"]["whole_catalog"] = dict(anchored_concordant=426, genes=28233, frac=426 / 28233,
                                             note="background includes 13,551 untested genes")
ue3 = ueval.merge(ce[["gene", "ave_control"]], on="gene", how="left")
assert ue3.g.tolist() == ueval.g.tolist() and ue3.ave_control.notna().all()
b3 = make_bins(ue3, "ave_control", 10)
D["D3"] = {s: strat_test(b3, gen_e, inE[s]) for s in LISTS}
joint = bins10 * 2 + deg_e.astype(int)
D["D4"] = {s: strat_test(joint, gen_e, inE[s]) for s in LISTS}
D5 = {}
for s in LISTS:
    t = TAB[(TAB.list == f"Kumak_{s}") & (TAB.match_state == "matched")]
    notdeg = t[~t.bulk_canonical_deg.astype(bool)]
    D5[s] = dict(matched_not_canonical_deg=len(notdeg), of_which_in_bulk_table=int(notdeg.in_bulk_table.sum()),
                 not_in_bulk_table=int((~notdeg.in_bulk_table.astype(bool)).sum()),
                 classes=notdeg.primary_evidence_class.value_counts().to_dict())
    if args.debug_null:
        D5[s]["note"] = "counts and classes from the real list; stratified test from the debug null list"
    lst5 = inE[s] & ~deg_e
    D5[s]["U_eval_test"] = strat_test(bins10, gen_e, lst5)
D["D5"] = D5
# D6: one genetic gene per 500 kb cluster
pp_ = pos.sort_values(["cluster", "susie", "g"], ascending=[True, False, True])
keep6 = set(pp_.drop_duplicates("cluster").g)
drop6 = set(pos.g) - keep6
u6 = ueval[~ueval.g.isin(drop6)].reset_index(drop=True)
b6 = make_bins(u6, "AveExpr", 10)
D6 = dict(clusters=int(pos.cluster.nunique()), U_eval_after_pruning=len(u6), genetic_after_pruning=int(u6.gen.sum()))
for s in LISTS:
    lg = set(pos.g[pos.g.isin(LSET[s])])
    D6[f"{s}_clusters_with_list_gene"] = int(pos.cluster[pos.g.isin(lg)].nunique())
    D6[f"{s}_test"] = strat_test(b6, u6.gen.to_numpy(), u6.g.isin(LSET[s]).to_numpy())
D["D6"] = D6
pos.assign(kept=pos.g.isin(keep6), in_NAS=pos.g.isin(LSET["NAS"]), in_Fibrosis=pos.g.isin(LSET["Fibrosis"]))[
    ["g", "susie", "susie_driving_gwas", "top_snp", "cluster", "kept", "in_NAS", "in_Fibrosis"]].to_csv(
    f"{OUT}/d6_clusters.tsv", sep="\t", index=False)
D7 = TAB[TAB.list_symbol.isin(NAMED) & TAB.list.str.startswith("Kumak")]
D7.to_csv(f"{OUT}/d7_named_genes.tsv", sep="\t", index=False)


def unmatched_reason(sym):
    if re.fullmatch(r"C\d+orf\d+", sym) or re.fullmatch(r"KIAA\d+", sym):
        return "retired alias"
    if re.fullmatch(r"[A-Z0-9]+-(AS|DT|IT|OT)\d*", sym):
        return "other"
    if re.fullmatch(r"[A-Z0-9]+-[A-Z0-9]+", sym):
        return "readthrough"
    if re.search(r"\dP\d*$", sym):
        return "pseudogene"
    return "other"


D8 = TAB[TAB.match_state == "unmatched"][["list", "list_symbol"]].copy()
D8["reason_class"] = D8.list_symbol.map(unmatched_reason)
D8.to_csv(f"{OUT}/d8_unmatched.tsv", sep="\t", index=False)
D["D8"] = {lst: grp.reason_class.value_counts().to_dict() for lst, grp in D8.groupby("list")}
D9 = TAB[TAB.list == "Elison"].copy()
ba = D9[D9.elison_finding_type == "boundary_absent"]
D["D9_boundary_absent"] = {
    r.list_symbol: (f"class {r.primary_evidence_class}; bulk {r.bulk_call_state} (log2FC {r.bulk_log2fc:.3f}, q {r.bulk_q:.3g}; "
                    f"canonical rule needs q < 0.05 and |log2FC| > 0.5); genetics {r.genetic_call_state} "
                    f"({r.genetic_testability_reason}; susie_state_t12 {r.susie_state_t12}, best SuSiE PP.H4 {float(r.coloc_best_susie_pp4):.3g})")
    for r in ba.itertuples()}
RES["descriptives"] = D

# ======================================================================= reading (prespec "Reading")
READ = {}
for s in LISTS:
    if LS[s]["untestable"]:
        READ[s] = dict(reading="untestable", rule=LS[s]["untestable"])
        continue
    r = T1[s]
    side = np.sign(r["O"] - r["E_strat"])
    d3s = bool(d3_ok and np.sign(D["D3"][s]["O"] - D["D3"][s]["E_strat"]) == side and side != 0)
    d6s = bool(np.sign(D6[f"{s}_test"]["O"] - D6[f"{s}_test"]["E_strat"]) == side and side != 0)
    conds = {"T1 Holm-adjusted p < 0.05": r["p_holm"] < 0.05, "PC-A": RES["PC_A"]["pass"],
             "PC-B": LS[s]["pcb_pass"], "PC-C": LS[s]["pcc_pass"], "NC-1 (S6)": not LS[s]["S6"],
             "D3 direction": d3s, "D6 direction": d6s}
    if all(conds.values()):
        rd = "directional_more" if side > 0 else "directional_less"
        failing = []
    elif conds["T1 Holm-adjusted p < 0.05"]:
        rd = "holm_significant_check_failed"
        failing = [k for k, v in conds.items() if not v]
    else:   # null reading: only PC-C and S6 qualify it (PC-A and PC-B already passed to get here)
        failing = [k for k in ("PC-C", "NC-1 (S6)") if not conds[k]]
        rd = "null" if not failing else "null_check_named"
    READ[s] = dict(reading=rd, conditions=conds, failing=failing)
RES["reading"] = READ
LONG = {"NAS": "NAS", "Fibrosis": "fibrosis"}
PLAIN = {"PC-C": "the positive control for the stratified test (PC-C) failed",
         "NC-1 (S6)": "the expression-matched random lists (NC-1) fell outside their limits",
         "D3 direction": "the control-only expression stratifier (D3) gave the other side of 1",
         "D6 direction": "one gene per signal cluster (D6) gave the other side of 1",
         "PC-A": "the reproduction check (PC-A) failed", "PC-B": "the planted-signal check (PC-B) failed"}


def plain_fail(s):
    return " and ".join(PLAIN[k] for k in READ[s]["failing"])


def fmt(r):
    s_ = f"O {r['O']}, E_strat {r['E_strat']:.3f}, O/E {r['O_over_E']:.3f}"
    if "MH_OR" in r:
        s_ += f", MH OR {r['MH_OR']:.3f} (95% CI {r['MH_lo']:.3f} to {r['MH_hi']:.3f})"
    return s_ + f", exact p {r['p_exact']:.4g}"


def list_verdict(s):
    rd = READ[s]["reading"]
    if rd == "untestable":
        d2 = D["D2"][s]
        return (f"{s}: untestable ({READ[s]['rule']}); expression not controlled: unstratified O {d2['O']}, "
                f"E0 {d2['E0']:.2f}, Fisher OR {d2['fisher_OR']:.3f}, nominal p {d2['fisher_p']:.4g}.")
    r = T1[s]
    core = (f"O {r['O']}, E_strat {r['E_strat']:.3f} ({r['K']} AveExpr bins), O/E_strat {r['O_over_E']:.3f}, "
            f"MH OR {r['MH_OR']:.3f} (RBG 95% CI {r['MH_lo']:.3f} to {r['MH_hi']:.3f}), exact two-sided T1 p {r['p_exact']:.4g}, "
            f"Holm-adjusted p {r['p_holm']:.4g} (family F1, 2 tests)")
    pcb = f"PC-B detection {LS[s]['pcb_up']['rate']:.2f} (OR up shift {SHIFT[s][0]:+d}) and {LS[s]['pcb_down']['rate']:.2f} (shift {SHIFT[s][1]:+d})"
    if rd.startswith("directional"):
        return f"{s}: colocalized {'more' if rd.endswith('more') else 'less'} often than genes of the same expression level. {core}."
    if rd == "holm_significant_check_failed":
        return f"{s}: T1 Holm-significant but no directional wording, because {plain_fail(s)}. {core}."
    if rd == "null":
        return f"{s}: null, the colocalization rate did not differ detectably from that of genes of the same expression level. {core}; {pcb}."
    return (f"{s}: null, did not differ detectably from the stratified expectation; {plain_fail(s)}, so the phrase "
            f"'of the same expression level' is not allowed. {core}; {pcb}.")


RES["lists"] = {s: {k: v for k, v in LS[s].items() if k not in ("bins", "n", "Gk", "Lk", "pmf")} for s in LISTS}
verdict = " ".join(list_verdict(s) for s in LISTS)
if args.debug_null:
    verdict = "DEBUG NULL (lists replaced; not real). " + verdict


# ======================================================================= proposed text (returned as data)
def disc_clause():
    parts = []
    for s in LISTS:
        rd = READ[s]["reading"]
        if rd == "untestable":
            parts.append((s, "could not be compared with genes of the same expression level", None))
            continue
        r = T1[s]
        ci = f"Mantel-Haenszel OR {r['MH_OR']:.2f}, 95% CI {r['MH_lo']:.2f} to {r['MH_hi']:.2f}"
        if rd == "directional_more":
            parts.append((s, "was higher than that of genes of the same expression level", ci))
        elif rd == "directional_less":
            parts.append((s, "was lower than that of genes of the same expression level", ci))
        elif rd == "null":
            parts.append((s, "did not differ detectably from that of genes of the same expression level", ci))
        elif rd == "null_check_named":
            parts.append((s, f"did not differ detectably from its stratified expectation, a comparison that does not fully "
                             f"control expression level because {plain_fail(s)}", ci))
        else:
            parts.append((s, f"differed from its stratified expectation, but {plain_fail(s)}, so no direction is stated", ci))
    if parts[0][1] == parts[1][1] and parts[0][2] and parts[1][2]:
        return (f"the colocalization rate {parts[0][1]} (NAS {parts[0][2]}; fibrosis "
                f"{parts[1][2].replace('Mantel-Haenszel ', '')})")
    out = []
    for s, txt, ci in parts:
        out.append(f"that of {LONG[s]} genes {txt}" + (f" ({ci})" if ci else ""))
    return "the colocalization rate of " + out[0][len("that of "):] + ", and " + out[1]


def detail_sentence():
    bits = []
    for s in LISTS:
        if READ[s]["reading"] == "untestable":
            d2 = D["D2"][s]
            bits.append(f"{LONG[s]}: {d2['O']} colocalized genes against {d2['E0']:.1f} expected without expression stratification "
                        f"(expression not controlled)")
            continue
        r = T1[s]
        bits.append(f"{LONG[s]}: {r['O']} colocalized genes against {r['E_strat']:.1f} expected; planted shifts equal to an OR of about "
                    f"{G3[s]['OR_up']:.2f} and {G3[s]['OR_dn']:.2f} were detected in {100 * LS[s]['pcb_up']['rate']:.0f}% and "
                    f"{100 * LS[s]['pcb_down']['rate']:.0f}% of simulations")
    return "In the evaluated genes, " + "; ".join(bits) + "."


R_N = D1["Kumak_NAS"]["supported"]; R_F = D1["Kumak_Fibrosis"]["supported"]
old_disc = "results are credited source findings rather than results of this paper."
new_disc = (old_disc + " Of the 944 NAS-associated and 1,013 fibrosis-associated genes they report that map to the Catalog, "
            f"{R_N:,} and {R_F:,} reach the canonical bulk call here; of those present in both the bulk and colocalization tables, "
            f"661 of 902 and 714 of 943 had no multi-signal colocalization result, and among the 241 and 228 evaluated genes "
            f"{disc_clause()}. {detail_sentence()}")
disc_lines = open(DISC).read().splitlines()
disc_ln = [i + 1 for i, l in enumerate(disc_lines) if "results are credited source findings rather than results of this" in l]
comp_lines = open(COMP).read().splitlines()
comp_ln = [i + 1 for i, l in enumerate(comp_lines) if "which shows the inherited axis is genuinely orthogonal" in l]
old_comp = ("we are not weaker there; and how little of their list touches the genetic map,\n"
            "   which shows the inherited axis is genuinely orthogonal and unreachable from\n   bulk expression alone.")
RPLAIN = {"directional_more": "colocalized more often than genes of the same expression level",
          "directional_less": "colocalized less often than genes of the same expression level",
          "null": "did not differ detectably from genes of the same expression level",
          "null_check_named": "did not differ detectably from the stratified expectation (expression not fully controlled)",
          "holm_significant_check_failed": "Holm-significant, no direction stated because a check failed"}
t1txt = "; ".join(f"{LONG[s]} genes {RPLAIN[READ[s]['reading']]}: O {T1[s]['O']} vs E_strat {T1[s]['E_strat']:.1f}, "
                  f"MH OR {T1[s]['MH_OR']:.2f} ({T1[s]['MH_lo']:.2f} to {T1[s]['MH_hi']:.2f}), Holm p {T1[s]['p_holm']:.3g}"
                  if s in T1 else f"{LONG[s]} genes untestable" for s in LISTS)
new_comp = ("we are not weaker there (H71: " + f"{R_N} of 944 NAS and {R_F} of 1,013 fibrosis genes bulk-supported); "
            "and how much of their list touches the genetic map (H71, COLOC r2: 661 of 902 NAS and 714 of 943 fibrosis genes "
            f"in U_full have no multi-signal COLOC result; among evaluated genes, {t1txt}). The 'orthogonal' framing is not supported.")
diag("11_text", "old_text found verbatim in Discussion / competitor note",
     f"{old_disc in open(DISC).read()} / {old_comp in open(COMP).read()}")
PROPOSED = [dict(file=DISC, line=disc_ln, old_text=old_disc, new_text=new_disc),
            dict(file=COMP, line=comp_ln, old_text=old_comp, new_text=new_comp)]
RES["proposed_text"] = PROPOSED
with open(f"{OUT}/proposed_text_changes.json", "w") as f:
    json.dump(PROPOSED, f, indent=1)

# ======================================================================= RESULT_H71.md
L_ = ["## Part A (annotation, no inference; supplementary table `partA_requalification.tsv`, summary `partA_summary.tsv`)", ""]
for lst, d in D1.items():
    L_.append(f"- {lst}: {d['matched']} matched; bulk supported (D1) {d['supported']} ({100 * d['frac']:.1f}%), by class "
              f"{jdump(d['by_class'])}, by Kumak direction {jdump(d['by_direction'])}; q < 0.05 with |log2FC| <= 0.5: {d['q05_small_effect']} "
              f"(by class {jdump(d['q05_small_effect_by_class'])}).")
sm = PARTA["summary"]
for lst in ["Kumak_NAS", "Kumak_Fibrosis", "Elison"]:
    t = sm[sm.list == lst]
    for q in ["primary class", "susie_state_t12 (all matched)", "susie_state_t12 (in U_full)", "anchored+concordant by susie_state_t12"]:
        tt = t[t.quantity == q]
        L_.append(f"- {lst} {q}: " + ", ".join(f"{a} {float(b):g}" for a, b in zip(tt.level, tt.value)))
L_ += ["", "## Part B, primary T1 (U_eval, 3,769 genes; family F1 = T1_NAS and T1_Fibrosis, Holm, family-wise alpha 0.05)", ""]
for s in LISTS:
    if s not in T1:
        L_.append(f"- {s}: untestable ({LS[s]['untestable']}).")
        continue
    r = T1[s]
    L_.append(f"- {s}: {fmt(r)}; Holm-adjusted p {r['p_holm']:.4g}; permutation check p {r['p_perm']:.4g} "
              f"(band {r['perm_band']:.4g}), permutation mean - E_strat {r['perm_mean_minus_E']:.4f}; S7 not triggered. "
              f"RBG CI excludes 1: {r['rbg_excludes_1']}. Not tested: Breslow-Day {r['BD_chi2']:.3f} (nominal p {r['BD_p']:.4g}), "
              f"CMH {r['CMH_chi2']:.3f} (nominal p {r['CMH_p']:.4g}).")
L_ += ["", "## Controls", "", f"- PC-A reproduction: pass ({len([d for d in DIAG if d[0] == '3_PC-A'])} printed checks)."]
for s in LISTS:
    if "pcb_up" in LS[s]:
        L_.append(f"- PC-B {s}: detection at exact T1 p < 0.025: up {SHIFT[s][0]:+d} {LS[s]['pcb_up']['rate']:.3f}, down {SHIFT[s][1]:+d} "
                  f"{LS[s]['pcb_down']['rate']:.3f} (predicted 0.60 to 0.85; S5 needs >= 0.50 each; failed swaps "
                  f"{LS[s]['pcb_up']['failed']} / {LS[s]['pcb_down']['failed']}).")
    L_.append(f"- PC-C {s} (U_full, list vs canonical DEG): {fmt(LS[s]['pcc'])}; pass {LS[s]['pcc_pass']} (predicted O/E > 2.0, p < 1e-10).")
    if "nc1" in LS[s]:
        n1 = LS[s]["nc1"]
        L_.append(f"- NC-1 {s}: mean stratified O/E {n1['mean_strat_OE']:.3f}, fraction p < 0.05 {n1['fpr']:.3f}, mean unstratified O/E "
                  f"{n1['mean_unstrat_OE']:.3f}, distinct O values {n1['n_distinct_O']} (range {n1['O_min']} to {n1['O_max']}), "
                  f"median largest gap {n1['max_gap_median']:.4f}, max |SMD| {n1['smd_maxabs']:.4f}; prediction met {LS[s]['nc1_pred_met']}; "
                  f"S6 triggered {LS[s]['S6']}.")
    n2 = LS[s]["nc2"]
    L_.append(f"- NC-2 {s}: observed agreement {n2['observed']:.4f}, permuted mean {n2['perm_mean']:.4f} (expected {n2['expected']:.4f}), "
              f"permuted max {n2['perm_max']:.4f}; prediction met {n2['pred_met']}; observed above every permuted value {n2['observed_above_all']}.")
    if "identity_ok" in LS[s]:
        L_.append(f"- {s} identity check {LS[s]['identity_ok']}, role-swap check {LS[s]['roleswap_ok']}; G2 bins meeting rule {LS[s]['ok10']} of 10.")
L_ += ["", "## Reading conditions", ""]
for s in LISTS:
    L_.append(f"- {s}: reading {READ[s]['reading']}" + (f"; conditions {jdump(READ[s].get('conditions', {}))}" if "conditions" in READ[s] else ""))
L_ += ["", "## Descriptives (nominal p values, outside family F1, never called significant)", ""]
for s in LISTS:
    d2 = D["D2"][s]
    L_.append(f"- D2 {s} unstratified U_eval: O {d2['O']}, E0 {d2['E0']:.2f}, Fisher OR {d2['fisher_OR']:.3f}, p {d2['fisher_p']:.4g}.")
cw = D["D2_catalog_wide"]
L_.append(f"- D2 Catalog-wide: anchored+concordant NAS {cw['NAS']['anchored_concordant']}/{cw['NAS']['matched']} ({100 * cw['NAS']['frac']:.1f}%), "
          f"Fibrosis {cw['Fibrosis']['anchored_concordant']}/{cw['Fibrosis']['matched']} ({100 * cw['Fibrosis']['frac']:.1f}%), whole Catalog "
          f"426/28,233 ({100 * 426 / 28233:.1f}%); this background includes 13,551 untested genes.")
for s in LISTS:
    L_.append(f"- D3 {s} control-only stratifier (138 control units, 10 bins): {fmt(D['D3'][s])}.")
for s in LISTS:
    L_.append(f"- D4 {s} AveExpr decile x canonical DEG (20 strata): {fmt(D['D4'][s])}.")
for s in LISTS:
    d5 = D5[s]
    L_.append(f"- D5 {s} list genes that are not canonical DEGs: {d5['matched_not_canonical_deg']} of matched "
              f"({d5['of_which_in_bulk_table']} in the bulk table, {d5['not_in_bulk_table']} not); classes {jdump(d5['classes'])}; "
              f"U_eval stratified: {fmt(d5['U_eval_test'])}.")
L_.append(f"- D6 clusters among {len(pos)} genetic genes at 500 kb: {D6['clusters']}; U_eval after pruning {D6['U_eval_after_pruning']} "
          f"genes, {D6['genetic_after_pruning']} genetic.")
for s in LISTS:
    L_.append(f"- D6 {s}: clusters with a genetic list gene {D6[f'{s}_clusters_with_list_gene']}; one gene per cluster: {fmt(D6[f'{s}_test'])}.")
for r in D7.itertuples():
    L_.append(f"- D7 {r.list} {r.list_symbol}: {r.kumak_direction}, {r.primary_evidence_class}, log2FC {r.bulk_log2fc:.3f}, q {r.bulk_q:.3g}, "
              f"bulk {r.bulk_call_state}, susie_state_t12 {r.susie_state_t12}, genetics {r.genetic_call_state}.")
L_.append("- D8 unmatched symbols by reason class: " + "; ".join(f"{a} {jdump(b)}" for a, b in D["D8"].items()) + " (file d8_unmatched.tsv).")
for r in D9.itertuples():
    L_.append(f"- D9 Elison {r.list_symbol} ({r.elison_finding_type}): {r.primary_evidence_class}; susie_state_t12 {r.susie_state_t12}.")
for g_, txt in D["D9_boundary_absent"].items():
    L_.append(f"- D9 boundary_absent {g_}: {txt}.")
L_ += ["", "## Proposed text changes (returned as data, not installed)", ""]
for p in PROPOSED:
    L_ += [f"- File `{p['file']}`, line(s) {p['line']}.", f"  - old_text: {p['old_text']!r}", f"  - new_text: {p['new_text']!r}"]
L_ += ["", "Dependence: every genetic number rests on COLOC r2 (susie_coloc_r2_20260924T053244Z) and must be rebuilt on the adopted "
       "COLOC rerun. Genes are the unit; they are not replicates, and gene-level p values are likely too narrow.",
       "", "Files: diagnostics.tsv, results.json, partA_requalification.tsv, partA_summary.tsv, bins_aveexpr10_Ueval.tsv, "
       "pcb_replicates.tsv, nc1_lists.tsv, nc2_perms.tsv, t1_perm_sims.tsv, d6_clusters.tsv, d7_named_genes.tsv, d8_unmatched.tsv, "
       "proposed_text_changes.json, control_expression.tsv, sessionInfo.txt, pip_freeze.txt, script_sha256.tsv."]
finish("completed", verdict, L_)
