"""DIG1 job 59. Replication of the GSE130970 ballooning link (section 38) in two microarray cohorts with
component-level histology. Written before any expression value of these cohorts was read.
Modality: log2 microarray intensities (Illumina HT-12 GPL14951; Affymetrix HuGene 1.1 ST GPL11532). A linear
model on within-cohort z-scores of log2 intensity is native to this modality; no count model is used.
Cohorts and histology:
  GSE89632 (Arendt 2015; HC, SS, NASH): steatosis %, lobular inflammation, ballooning, fibrosis recorded.
  GSE48452 (Ahrens 2013): fat %, inflammation, NAS, fibrosis. Ballooning = NAS - steatosis grade - inflammation.
    The derivation is accepted only if every derived value is an integer in 0..2; otherwise untestable.
    After-surgery biopsies are excluded so that each patient contributes one biopsy.
  Steatosis grade from %: <=5 -> 0, >5-33 -> 1, >33-66 -> 2, >66 -> 3 (GSE48452 records 5% with NAS 0).
  Fibrosis stage used as recorded (GSE48452 half stages kept as numeric).
Expression: per symbol, mean log2 value over probes (log2 applied only if the matrix max > 100); within-cohort z
per gene; score = mean z over genes present. Sets as job 58: INFLAMED_HEP (FCAMR, TRIM31, CXCL10, IL32, CCL20),
IMMUNE (PTPRC, CD3E, CD2, CD68, CD163), HEPMASS (ALB, APOA1, CYP2E1, HP).
Model: score ~ steatosis + lobular + ballooning + fibrosis [+ IMMUNE + HEPMASS], complete cases, OLS, HC3 errors,
t on residual df.
Primary replication family: the ballooning coefficient of INFLAMED_HEP in the adjusted model, one test per
testable cohort, two-sided, Holm over the testable cohorts. Replication = positive and Holm p < 0.05.
Descriptive: unadjusted model; IMMUNE outcome; CCL20 and CXCL10 alone; fixed-effect inverse-variance summary of
the adjusted ballooning coefficient with GSE130970 (+0.238, SE 0.083, job 58b).
Unit: one biopsy per patient. Cross-sectional.
Usage: python -I 59_microarray_ballooning_replication.py <download_dir> <new_out_dir>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

assert len(sys.argv) == 3
DL, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=False)
SETS = {"INFLAMED_HEP": ["FCAMR", "TRIM31", "CXCL10", "IL32", "CCL20"], "IMMUNE": ["PTPRC", "CD3E", "CD2", "CD68", "CD163"],
        "HEPMASS": ["ALB", "APOA1", "CYP2E1", "HP"], "g_CCL20": ["CCL20"], "g_CXCL10": ["CXCL10"]}
COMP = ["steatosis", "lobular", "ballooning", "fibrosis"]


def read_matrix(f):
    meta, start = {}, None
    lines = gzip.open(f, "rt").read().split("\n")
    char = []
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            start = i + 1
            break
        f0 = line.split("\t")[0]
        v = [x.strip('"') for x in line.split("\t")[1:]]
        if f0 == "!Sample_characteristics_ch1":
            char.append(v)
        elif f0 in ("!Sample_geo_accession", "!Sample_title", "!Sample_source_name_ch1"):
            meta[f0] = v
    M = pd.DataFrame(meta)
    for v in char:
        for j, x in enumerate(v):
            if ": " in x:
                k, val = x.split(": ", 1)
                M.loc[j, k] = val
    end = next(i for i in range(start, len(lines)) if lines[i].startswith("!series_matrix_table_end"))
    rows = [l.split("\t") for l in lines[start:end]]
    hdr = [x.strip('"') for x in rows[0]]
    X = pd.DataFrame([[x.strip('"') for x in r] for r in rows[1:]], columns=hdr).set_index("ID_REF")
    X = X.apply(pd.to_numeric, errors="coerce")
    assert list(X.columns) == list(M["!Sample_geo_accession"])
    return M, X


def annot(f, kind):
    import io
    lines = open(f).read().split("\n")
    i0 = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")) + 1
    i1 = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_end"))
    A = pd.read_csv(io.StringIO("\n".join(lines[i0:i1])), sep="\t", dtype=str, low_memory=False)
    if kind == "illumina":
        return dict(zip(A.ID, A.Symbol))
    sym = A.gene_assignment.fillna("---").str.split(" // ").str[1].fillna("---").str.strip()
    return dict(zip(A.ID, sym))


def gene_matrix(X, p2s):
    if np.nanmax(X.values) > 100:
        X = np.log2(X.clip(lower=1))
    s = X.index.map(lambda p: p2s.get(p, None))
    X = X[pd.notna(s) & (s != "---")]
    return X.groupby(s[pd.notna(s) & (s != "---")]).mean()


def scores(G):
    Z = G.sub(G.mean(1), axis=0).div(G.std(1, ddof=1), axis=0)
    out, found = {}, {}
    for k, gs in SETS.items():
        present = [g for g in gs if g in Z.index]
        found[k] = present
        out[k] = Z.loc[present].mean(0) if present else pd.Series(np.nan, index=Z.columns)
    return pd.DataFrame(out), found


def hc3(y, Xd, names):
    X = np.column_stack([np.ones(len(y)), Xd]); n, k = X.shape
    XtXi = np.linalg.inv(X.T @ X); b = XtXi @ X.T @ y; e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    V = XtXi @ (X * (e / (1 - h))[:, None]).T @ (X * (e / (1 - h))[:, None]) @ XtXi
    se = np.sqrt(np.diag(V)); t = b / se
    return pd.DataFrame(dict(term=["(Intercept)"] + names, estimate=b, se=se, p=2 * stats.t.sf(np.abs(t), n - k), n=n))


def steat_grade(pct):
    pct = pd.to_numeric(pct, errors="coerce")
    return pd.cut(pct, [-np.inf, 5, 33, 66, np.inf], labels=[0, 1, 2, 3]).astype(float)


rows, status, found_all, dfs = [], [], {}, {}
# GSE89632
M, X = read_matrix(DL / "GSE89632_series_matrix.txt.gz")
G = gene_matrix(X, annot(DL / "GPL14951_annot.txt", "illumina"))
S, found_all["GSE89632"] = scores(G)
H = pd.DataFrame({"gsm": M["!Sample_geo_accession"], "dx": M["diagnosis"],
                  "steatosis": steat_grade(M["steatosis (%)"]),
                  "lobular": pd.to_numeric(M["lobular inflammation (severity)"], errors="coerce"),
                  "ballooning": pd.to_numeric(M["ballooning (intensity)"], errors="coerce"),
                  "fibrosis": pd.to_numeric(M["fibrosis (stage)"], errors="coerce"),
                  "nas": pd.to_numeric(M["nafld activity score"], errors="coerce")})
H["nas_check"] = H.steatosis + H.lobular + H.ballooning - H.nas
dfs["GSE89632"] = H.join(S, on="gsm")
status.append(dict(cohort="GSE89632", testable=True, note=f"NAS recomputed from components minus recorded NAS: "
                   f"{H.nas_check.value_counts(dropna=True).to_dict()}"))
# GSE48452
M, X = read_matrix(DL / "GSE48452_series_matrix.txt.gz")
G = gene_matrix(X, annot(DL / "GPL11532_annot.txt", "affy"))
S, found_all["GSE48452"] = scores(G)
H = pd.DataFrame({"gsm": M["!Sample_geo_accession"], "dx": M["group"], "surgery": M["bariatric surgery"],
                  "steatosis": steat_grade(M["fat"]), "lobular": pd.to_numeric(M["inflammation"], errors="coerce"),
                  "fibrosis": pd.to_numeric(M["fibrosis"], errors="coerce"), "nas": pd.to_numeric(M["nas"], errors="coerce")})
H = H[H.surgery != "after surgery"].copy()
H["ballooning"] = H.nas - H.steatosis - H.lobular
d = H.ballooning.dropna()
ok = bool(((d >= 0) & (d <= 2) & (d == d.round())).all())
status.append(dict(cohort="GSE48452", testable=ok, note=f"derived ballooning values {d.value_counts().to_dict()}; "
                   f"after-surgery excluded; n={len(H)}"))
dfs["GSE48452"] = H.join(S, on="gsm")

for co, D in dfs.items():
    D.to_csv(OUT / f"sample_table_{co}.tsv", sep="\t", index=False)
    st = next(s for s in status if s["cohort"] == co)
    if not st["testable"]:
        continue
    for outcome in SETS:
        for adj in (False, True):
            if outcome in ("IMMUNE", "HEPMASS") and adj:
                continue
            cols = COMP + (["IMMUNE", "HEPMASS"] if adj else [])
            dd = D.dropna(subset=cols + [outcome])
            r = hc3(dd[outcome].to_numpy(float), dd[cols].to_numpy(float), cols)
            r["cohort"], r["outcome"], r["adjusted"] = co, outcome, adj
            rows.append(r)
R = pd.concat(rows, ignore_index=True)
R = R[R.term != "(Intercept)"]
R.to_csv(OUT / "component_models.tsv", sep="\t", index=False)
P = R[(R.outcome == "INFLAMED_HEP") & R.adjusted & (R.term == "ballooning")].copy()
m = len(P)
order = np.argsort(P.p.to_numpy())
holm = np.empty(m); run = 0
for rank, i in enumerate(order):
    run = max(run, min(1, (m - rank) * P.p.to_numpy()[i])); holm[i] = run
P["p_holm"] = holm
P["replicates"] = (P.estimate > 0) & (P.p_holm < 0.05)
P.to_csv(OUT / "primary_replication.tsv", sep="\t", index=False)
E = pd.concat([P[["cohort", "estimate", "se"]], pd.DataFrame([dict(cohort="GSE130970", estimate=0.23775, se=0.0831)])])
w = 1 / E.se ** 2; est = (w * E.estimate).sum() / w.sum(); se = np.sqrt(1 / w.sum()); Q = (w * (E.estimate - est) ** 2).sum()
FE = pd.DataFrame([dict(k=len(E), estimate=est, se=se, p=2 * stats.norm.sf(abs(est / se)), Q=Q,
                        I2=max(0, (Q - (len(E) - 1)) / Q) if Q > 0 else 0)])
FE.to_csv(OUT / "fixed_effect_ballooning.tsv", sep="\t", index=False)
pd.DataFrame(status).to_csv(OUT / "cohort_status.tsv", sep="\t", index=False)
with open(OUT / "set_genes_found.txt", "w") as fh:
    for co, f in found_all.items():
        for k, v in f.items():
            fh.write(f"{co}\t{k}\t{','.join(v)}\n")
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250)
print(pd.DataFrame(status).to_string()); print(found_all)
for co, D in dfs.items():
    print(co, "ballooning x dx"); print(pd.crosstab(D.ballooning, D.dx))
print(P.to_string()); print(FE.to_string())
print(R[R.term.isin(COMP)].round(4).to_string())
