#!/usr/bin/env python3
"""
72_ag_multimodal_mechanism.py  --  Phase-4 AlphaGenome MULTIMODAL MECHANISM profile.

Uses the FULL AlphaGenome head set (not just the one RNA-direction channel that
src/65 called) to produce, per variant, a mechanism-CLASS / which-gene / magnitude
profile. It is the crown-jewel upgrade of the seqfunc which-gene axis: a 3D
contact-map linker that replaces the barely-orthogonal nearest-TSS baseline for the
convergent nominations.

STRICT GUARDRAILS (non-negotiable, enforced by construction):
  * Mechanism-CLASS / which-gene / magnitude ONLY. NEVER a direction claim -- the
    eQTL-direction gate FAILED for BOTH sequence models (Borzoi auROC 0.559 /
    AlphaGenome 0.561, ns). Signed quantiles are RECORDED but never turned into an
    up/down assertion.
  * NEVER a scored convergence channel (circularity vs COLOC + epigenomic). Outputs
    are absent from convergence_evidence.csv inputs; src/70 folds them in as Axis-5
    class/which-gene ANNOTATIONS only.
  * Contact linking = a which-gene / regulatory-loop CLASS hypothesis (2 kb-binned,
    tissue-AVERAGED contact map -- only 1/28 AG contact tracks is liver/HepG2, and
    TAD/loop architecture is largely cell-type-invariant; the variant-delta is
    unbenchmarked and used class-only). Never a quantitative contact-strength claim.
  * Corroborators (ABC, Roadmap E066, SCREEN cCRE, ChromBPNet, motifbreakR) are
    MEASURED + INDEPENDENT -> reported as agreement FLAGS, never summed into a score.
    AG-contact is scored against Borzoi + ABC (independent model/modality), NEVER
    against AG's own expression head.
  * AG outputs must not train other models (license); research-only.

FIVE COMPONENTS
  C1 Fingerprint : one score_variant call/variant with the 11 scalar scorers
     (RNA_SEQ, POLYADENYLATION, SPLICE_JUNCTIONS/SITES/SITE_USAGE, ATAC, DNASE,
     CHIP_HISTONE, CHIP_TF, CAGE, PROCAP), liver-masked; per head -> signed-extreme
     |quantile|. histone_mark / transcription_factor tidy cols -> element identity.
  C2 Contact which-gene : predict_variant([CONTACT_MAPS]) -> REF+ALT 512x512x28 @
     2048 bp. Tissue-average REF map; for each in-window candidate TSS read
     contact(variant_bin, TSS_bin); ag_contact_gene = argmax; contact_rewire_abs =
     |ALT-REF| at the winning cell (class-only). HepG2-only argmax as an agreement
     flag.
  C3 Mechanism class : ordered rule over the comparable liver |quantiles| ->
     primary + secondary class {splice-altering, promoter/TSS, enhancer-disrupting,
     TF-footprint-breaking, 3D-contact-rewiring} + confidence (argmax margin).
  C4 Corroboration : ag_contact_gene vs ABC + Borzoi; element class vs Roadmap E066
     + SCREEN cCRE; TF vs ChromBPNet + motifbreakR. Flags only.

ADDITIVE ONLY -- new file (script 72). Reuses src/65 liver-mask + key/retry/backoff
and src/64b GTF/MANE-TSS resolution (imported by path). Writes only under
results/seqfunc/. NO GPU (AlphaGenome is server-side); io/cpu node.

POLARITY CONTRACT (identical to 64b/65): score strictly hg38 REF -> hg38 ALT. The
regulatory leads and coding effectors already carry hg38-oriented alleles
(ref_hg38/alt_hg38, variant_hg38); FASTA is used only for an optional base check.
"""

import argparse
import csv
import importlib.util
import os
import re
import sys
import time

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ constants
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
LIVER_ONTOLOGY = {"UBERON:0002107", "CL:0000182"}
LIVER_NAME_RE = re.compile(r"\bliver\b|hepatocyte|hepg2|hepatic", re.IGNORECASE)

# scalar heads scored in ONE score_variant call (verified live: 11 keys OK)
SCALAR_SCORER_KEYS = ["RNA_SEQ", "POLYADENYLATION",
                      "SPLICE_JUNCTIONS", "SPLICE_SITES", "SPLICE_SITE_USAGE",
                      "ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "PROCAP"]

# hg38 anchors (verified against the reference FASTA in 64b/65).
ANCHORS = [
    {"variant_id": "SORT1_rs12740374", "chr": "1", "pos_hg38": 109274968,
     "hg38_ref": "G", "hg38_alt": "T", "gene": "SORT1",
     "ensembl": "ENSG00000134243", "anchor_name": "SORT1_posctrl",
     "variant_id_hg19": "1:109817490:G:T"},
    {"variant_id": "HSD17B13_rs72613567", "chr": "4", "pos_hg38": 87310240,
     "hg38_ref": "T", "hg38_alt": "TA", "gene": "HSD17B13",
     "ensembl": "ENSG00000170509", "anchor_name": "HSD17B13_splice",
     "variant_id_hg19": "4:88231391:T:TA"},
    {"variant_id": "PNPLA3_rs738409", "chr": "22", "pos_hg38": 43928847,
     "hg38_ref": "C", "hg38_alt": "G", "gene": "PNPLA3",
     "ensembl": "ENSG00000100344", "anchor_name": "PNPLA3_I148M",
     "variant_id_hg19": "22:44324727:C:G"},
]

# Anchor competitor candidate TSS (so the contact which-gene test is discriminating,
# not a single-gene tautology). SORT1 1p13 locus: PSRC1 (9 kb) + CELSR2 (25 kb) are the
# NEAREST genes; SORT1 is the DISTAL (123 kb) causal hepatic effector (Musunuru 2010) --
# contact must overturn nearest-TSS to recover SORT1. Ensembl ids + presence verified in
# GENCODE v49.
ANCHOR_CANDIDATES = {
    "SORT1_posctrl": [("CELSR2", "ENSG00000143126"), ("PSRC1", "ENSG00000134222"),
                      ("SORT1", "ENSG00000134243"), ("PSMA5", "ENSG00000143106"),
                      ("SYPL2", "ENSG00000143028"), ("ATXN7L2", "ENSG00000162650"),
                      ("AMIGO1", "ENSG00000181754")],
    "HSD17B13_splice": [("HSD17B13", "ENSG00000170509"), ("MTTP", "ENSG00000138823"),
                        ("HSD17B11", "ENSG00000198189")],
    "PNPLA3_I148M": [("PNPLA3", "ENSG00000100344"), ("SAMM50", "ENSG00000100307"),
                     ("PARVB", "ENSG00000188677")],
}

# mechanism-class floor on the MEAN liver |quantile| (mean aggregation centres the
# null near 0; 0.50 = the class's typical liver signal is a top-half calibrated effect).
CLASS_ACTIVE_Q = 0.50
# a contact rewire is "high" if in the top cohort quintile (class-only, relative).
CONTACT_HIGH_PCTILE = 80.0
# promoter/enhancer split by variant-to-TSS distance (bp).
PROMOTER_TSS_BP = 2000

SF = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc"
GWAS_ATAC = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/gwas_atac"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def strip_ver(ens):
    return ens.split(".")[0] if ens else ens


# ------------------------------------------------------------------ import 64b TSS machinery
def load_64b(src_dir):
    path = os.path.join(src_dir, "64b_borzoi_eqtl_logsed.py")
    spec = importlib.util.spec_from_file_location("borzoi64b", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ key handling (from src/65)
def load_api_key(key_file):
    k = os.environ.get("ALPHA_GENOME_API_KEY")
    if k and k.strip():
        return k.strip(), "env:ALPHA_GENOME_API_KEY"
    kf = key_file or os.path.expanduser("~/.alphagenome_key")
    if os.path.exists(kf):
        with open(kf) as f:
            v = f.read().strip()
        if v:
            return v, f"keyfile:{kf}"
    return None, None


# ------------------------------------------------------------------ liver mask (from src/65)
def liver_mask(df):
    m = np.zeros(len(df), dtype=bool)
    if "ontology_curie" in df.columns:
        m |= df["ontology_curie"].astype(str).isin(LIVER_ONTOLOGY).values
    for col in ("gtex_tissue", "biosample_name", "track_name", "name"):
        if col in df.columns:
            m |= df[col].astype(str).str.contains(LIVER_NAME_RE).values
    return m


def signed_extreme(series):
    """Signed value with the largest |.| (NaN-safe). Magnitude used for class;
    sign RECORDED but never asserted as direction."""
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return float("nan")
    return float(s.iloc[int(np.abs(s.values).argmax())])


# ------------------------------------------------------------------ fingerprint extraction
def head_quantile(df, output_type, gene_ensembl=None, gene_symbol=None,
                  scorer_substr=None, scorer_exclude=None):
    """Liver-masked signed-extreme quantile for one head, optionally gene-filtered
    (RNA gene LFC) and scorer-filtered (split PolyA vs GeneMaskLFC inside RNA_SEQ)."""
    if "output_type" not in df.columns:
        return float("nan"), 0
    sub = df[df["output_type"].astype(str) == output_type]
    if scorer_substr is not None and "variant_scorer" in sub.columns:
        sub = sub[sub["variant_scorer"].astype(str).str.contains(scorer_substr, case=False, na=False)]
    if scorer_exclude is not None and "variant_scorer" in sub.columns:
        sub = sub[~sub["variant_scorer"].astype(str).str.contains(scorer_exclude, case=False, na=False)]
    if len(sub) == 0:
        return float("nan"), 0
    if gene_ensembl is not None:
        sid = strip_ver(gene_ensembl)
        gm = np.zeros(len(sub), dtype=bool)
        if "gene_id" in sub.columns:
            gm |= sub["gene_id"].astype(str).map(strip_ver).values == sid
        if gene_symbol is not None and "gene_name" in sub.columns:
            gm |= (sub["gene_name"].astype(str) == gene_symbol).values
        subg = sub[gm]
        sub = subg if len(subg) else sub
    lm = liver_mask(sub)
    subl = sub[lm] if lm.sum() > 0 else sub
    # MEAN signed quantile over liver tracks -- calibrated + multiplicity-robust.
    # (signed-EXTREME saturates: the max |quantile| over N calibrated tracks -> ~1.0
    # for any variant when N is large, e.g. 189 splice-junction or 560 TF tracks.)
    q = pd.to_numeric(subl.get("quantile_score", pd.Series(dtype=float)), errors="coerce").dropna()
    return (float(q.mean()) if len(q) else float("nan")), int(lm.sum())


def top_element(df, output_type, id_col):
    """Top liver track's identity (transcription_factor / histone_mark) by |quantile|."""
    if "output_type" not in df.columns or id_col not in df.columns:
        return "", float("nan")
    sub = df[df["output_type"].astype(str) == output_type].copy()
    if len(sub) == 0:
        return "", float("nan")
    lm = liver_mask(sub)
    sub = sub[lm] if lm.sum() > 0 else sub
    sub = sub[sub[id_col].astype(str).str.strip().replace({"nan": "", "None": ""}) != ""]
    if len(sub) == 0:
        return "", float("nan")
    sub = sub.assign(_aq=pd.to_numeric(sub["quantile_score"], errors="coerce").abs())
    sub = sub.dropna(subset=["_aq"]).sort_values("_aq", ascending=False)
    if len(sub) == 0:
        return "", float("nan")
    r = sub.iloc[0]
    return str(r[id_col]), float(pd.to_numeric(r["quantile_score"], errors="coerce"))


def histone_mark_quantile(df, mark_substr):
    """Mean liver signed quantile among CHIP_HISTONE rows of one histone_mark
    (per-mark so a promoter H3K4me3 signal is not diluted by unrelated marks;
    mean over that mark's liver tracks, multiplicity-robust)."""
    if "output_type" not in df.columns or "histone_mark" not in df.columns:
        return float("nan")
    sub = df[df["output_type"].astype(str) == "CHIP_HISTONE"]
    sub = sub[sub["histone_mark"].astype(str).str.contains(mark_substr, case=False, na=False)]
    if len(sub) == 0:
        return float("nan")
    lm = liver_mask(sub)
    sub = sub[lm] if lm.sum() > 0 else sub
    q = pd.to_numeric(sub.get("quantile_score", pd.Series(dtype=float)), errors="coerce").dropna()
    return float(q.mean()) if len(q) else float("nan")


def extract_fingerprint(df, target_ens, target_sym):
    """Flatten one variant's tidy_scores into the fingerprint fields."""
    fp = {}
    # gene expression LFC (GeneMaskLFCScorer) on the target gene, liver
    fp["q_rna_gene"], fp["n_liver_rna"] = head_quantile(
        df, "RNA_SEQ", gene_ensembl=target_ens, gene_symbol=target_sym,
        scorer_exclude="poly")
    fp["q_polya"], _ = head_quantile(df, "RNA_SEQ", gene_ensembl=target_ens,
                                     gene_symbol=target_sym, scorer_substr="poly")
    for key, ot in (("q_splice_junctions", "SPLICE_JUNCTIONS"),
                    ("q_splice_sites", "SPLICE_SITES"),
                    ("q_splice_site_usage", "SPLICE_SITE_USAGE"),
                    ("q_atac", "ATAC"), ("q_dnase", "DNASE"),
                    ("q_chip_histone", "CHIP_HISTONE"), ("q_chip_tf", "CHIP_TF"),
                    ("q_cage", "CAGE"), ("q_procap", "PROCAP")):
        fp[key], _ = head_quantile(df, ot)
    fp["top_tf"], fp["top_tf_quantile"] = top_element(df, "CHIP_TF", "transcription_factor")
    fp["top_histone_mark"], fp["top_histone_quantile"] = top_element(df, "CHIP_HISTONE", "histone_mark")
    fp["q_h3k4me3"] = histone_mark_quantile(df, "H3K4me3")
    fp["q_h3k27ac"] = histone_mark_quantile(df, "H3K27ac")
    fp["q_h3k4me1"] = histone_mark_quantile(df, "H3K4me1")
    return fp


# ------------------------------------------------------------------ mechanism class
def _absnan(x):
    try:
        v = abs(float(x))
        return v if v == v else 0.0
    except Exception:
        return 0.0


# Roadmap 15-state ChromHMM (core model) numeric codes -> element class.
# E1 TssA, E2 TssAFlnk, E10 TssBiv = promoter; E6 EnhG, E7 Enh, E12 EnhBiv = enhancer.
PROMOTER_STATES = {"E1", "E2", "E10"}
ENHANCER_STATES = {"E6", "E7", "E12"}
CHROMHMM_NAMES = {"E1": "TssA", "E2": "TssAFlnk", "E3": "TxFlnk", "E4": "Tx",
                  "E5": "TxWk", "E6": "EnhG", "E7": "Enh", "E8": "ZNF/Rpts",
                  "E9": "Het", "E10": "TssBiv", "E11": "BivFlnk", "E12": "EnhBiv",
                  "E13": "ReprPC", "E14": "ReprPCWk", "E15": "Quies"}


def assign_mechanism_class(fp, tss_dist, contact_high, tf_corroborated=False,
                           chromhmm_state=""):
    """Ordered rule over MEAN liver |quantiles| -> primary+secondary class + confidence.
    CLASS, never direction. Uses only discriminative heads (splice_junctions is dropped
    -- its null is degenerate, saturating for any variant). TF-footprint is asserted only
    when an independent footprint caller (ChromBPNet / motifbreakR) corroborates; the
    MEASURED Roadmap chromHMM state disambiguates promoter vs enhancer."""
    splice = max(_absnan(fp.get("q_splice_sites")), _absnan(fp.get("q_splice_site_usage")))
    prom = max(_absnan(fp.get("q_cage")), _absnan(fp.get("q_procap")), _absnan(fp.get("q_h3k4me3")))
    if tss_dist is not None and tss_dist > PROMOTER_TSS_BP:
        prom *= 0.5   # distal -> down-weight the promoter hypothesis
    enh = max(_absnan(fp.get("q_atac")), _absnan(fp.get("q_dnase")),
              _absnan(fp.get("q_h3k27ac")), _absnan(fp.get("q_h3k4me1")))
    tf = _absnan(fp.get("q_chip_tf"))
    scores = {"splice-altering": splice, "promoter/TSS": prom,
              "enhancer-disrupting": enh, "TF-footprint-breaking": tf}
    ordered = sorted(scores.items(), key=lambda kv: -kv[1])
    # TF-footprint requires an independent footprint vote; else drop it from contention
    if ordered[0][0] == "TF-footprint-breaking" and not tf_corroborated:
        ordered = [kv for kv in ordered if kv[0] != "TF-footprint-breaking"]
    primary, p_score = ordered[0]
    secondary, s_score = ordered[1]
    # MEASURED chromHMM promoter/enhancer override (more trustworthy than saturated heads)
    st = str(chromhmm_state).strip()
    if primary in ("promoter/TSS", "enhancer-disrupting"):
        if st in PROMOTER_STATES:
            primary = "promoter/TSS"
        elif st in ENHANCER_STATES:
            primary = "enhancer-disrupting"
    # 3D-contact-rewiring override: local heads quiet but a strong long-range loop
    if p_score < CLASS_ACTIVE_Q and contact_high:
        return ("3D-contact-rewiring", primary, round(p_score, 4), scores)
    if p_score < CLASS_ACTIVE_Q:
        primary = "sub-threshold"
    return (primary, secondary, round(p_score - s_score, 4), scores)


# ------------------------------------------------------------------ corroborator loaders
def load_abc():
    """hg19 variant_id -> set of ABC target genes (Nasser 2021 liver AvgHiC)."""
    p = os.path.join(GWAS_ATAC, "abc_variant_to_gene.csv")
    d = {}
    if not os.path.exists(p):
        return d
    for r in csv.DictReader(open(p)):
        v = str(r.get("variant_id", "")).replace("chr", "")
        g = str(r.get("abc_target_gene", "")).strip()
        if v and g:
            d.setdefault(v, set()).add(g)
    return d


def load_abc_bed(path):
    """hg19 ABC liver enhancer BED (chr,start,end,TargetGene,...) -> {chr: sorted
    [(start,end,gene)]}. Universe-INDEPENDENT ABC which-gene corroborator: overlap the
    lead's hg19 position against enhancer elements (Nasser 2021 liver/HepG2 AvgHiC).
    Complements the sparse precomputed abc_variant_to_gene.csv (built on 55b's universe)."""
    d = {}
    if not path or not os.path.exists(path):
        return d
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci, si, ei, gi = (header.index("chr"), header.index("start"),
                              header.index("end"), header.index("TargetGene"))
        except ValueError:
            ci, si, ei, gi = 0, 1, 2, 3
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= gi:
                continue
            ch = f[ci].replace("chr", "")
            try:
                d.setdefault(ch, []).append((int(f[si]), int(f[ei]), f[gi]))
            except ValueError:
                continue
    for ch in d:
        d[ch].sort()
    return d


def abc_bed_genes(abc_bed, chrom, pos1_hg19):
    """Set of ABC TargetGenes whose enhancer element spans the hg19 position."""
    recs = abc_bed.get(str(chrom).replace("chr", ""), [])
    out = set()
    for (s, e, g) in recs:
        if s <= pos1_hg19 <= e:
            out.add(g)
        elif s > pos1_hg19:
            break
    return out


def load_ccre():
    """hg19 variant_id -> (ccre_class_simple, gene_assigned)."""
    p = os.path.join(GWAS_ATAC, "ccre_variant_overlap.csv")
    d = {}
    if not os.path.exists(p):
        return d
    for r in csv.DictReader(open(p)):
        v = str(r.get("variant_id", "")).replace("chr", "")
        d[v] = (str(r.get("ccre_class_simple", r.get("ccre_class", ""))),
                str(r.get("gene_assigned", "")))
    return d


def load_chrombpnet():
    """hg19 variant_id -> disrupted_tf_motif (ChromBPNet HepG2 footprint)."""
    p = os.path.join(SF, "chrombpnet_accessibility.tsv")
    d = {}
    if not os.path.exists(p):
        return d
    for r in csv.DictReader(open(p), delimiter="\t"):
        v = str(r.get("variant_id_hg19", r.get("variant", ""))).replace("chr", "")
        tf = str(r.get("disrupted_tf_motif", "")).strip()
        if v and tf and tf != "-":
            d[v] = tf
    return d


def load_motifbreakr():
    """hg19 SNP_id -> set of disrupted TF names (motifbreakR + FIMO)."""
    p = os.path.join(GWAS_ATAC, "motif_disruption_scores.csv")
    d = {}
    if not os.path.exists(p):
        return d
    for r in csv.DictReader(open(p)):
        v = str(r.get("SNP_id", "")).replace("chr", "")
        tf = str(r.get("tf_name", "")).strip()
        if v and tf:
            d.setdefault(v, set()).add(tf)
    return d


def load_roadmap_beds(roadmap_dir):
    """chromHMM state / H3K27ac / H3K4me3 hg38 interval trees for liver E066."""
    beds = {}
    try:
        import gzip
        def read_bed(fn, name_col=3):
            path = os.path.join(roadmap_dir, fn)
            recs = {}
            if not os.path.exists(path):
                return recs
            op = gzip.open if path.endswith(".gz") else open
            with op(path, "rt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= name_col:
                        name = "peak"
                    else:
                        name = f[name_col]
                    ch = f[0]
                    recs.setdefault(ch, []).append((int(f[1]), int(f[2]), name))
            for ch in recs:
                recs[ch].sort()
            return recs
        beds["chromHMM"] = read_bed("E066_15_chromHMM_hg38.bed.gz", 3)
        beds["H3K27ac"] = read_bed("E066_H3K27ac_hg38.broadPeak.gz", 3)
        beds["H3K4me3"] = read_bed("E066_H3K4me3_hg38.narrowPeak.gz", 3)
    except Exception as e:  # noqa
        log(f"roadmap load warning: {e}")
    return beds


def roadmap_hit(beds, kind, chrom, pos1):
    ch = chrom if chrom.startswith("chr") else "chr" + chrom
    recs = beds.get(kind, {}).get(ch, [])
    for (s, e, name) in recs:
        if s <= pos1 <= e:
            return name
        if s > pos1:
            break
    return ""


# ------------------------------------------------------------------ candidate TSS resolution
def build_candidate_index(candidates_path):
    """variant_id_hg19 -> list of {gene, ensembl, biotype, tss_dist, is_best,
    magnitude_percentile}. Source: borzoi_magnitude_candidates.tsv."""
    idx = {}
    if not os.path.exists(candidates_path):
        return idx
    for r in csv.DictReader(open(candidates_path), delimiter="\t"):
        v = str(r.get("variant_id_hg19", "")).replace("chr", "")
        idx.setdefault(v, []).append({
            "gene": r.get("gene", ""), "ensembl": r.get("ensembl", ""),
            "biotype": r.get("biotype", ""),
            "tss_dist": r.get("tss_dist", ""), "is_best": r.get("is_best", ""),
            "magnitude_percentile": r.get("magnitude_percentile", ""),
        })
    return idx


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leads", default=os.path.join(SF, "borzoi_magnitude_leads.tsv"))
    ap.add_argument("--candidates", default=os.path.join(SF, "borzoi_magnitude_candidates.tsv"))
    ap.add_argument("--nominations", default=os.path.join(SF, "borzoi_magnitude_nominations.tsv"))
    ap.add_argument("--coding", default=os.path.join(SF, "coding_hardening.tsv"))
    ap.add_argument("--substrate", default=os.path.join(SF, "variant_substrate_hg38.tsv"))
    ap.add_argument("--gtf", default="/gpfs/commons/home/jameslee/reference_genome/"
                    "gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
    ap.add_argument("--fasta", default="/gpfs/commons/home/jameslee/reference_genome/"
                    "refdata-gex-GRCh38-2024-A/fasta/genome.fa")
    ap.add_argument("--roadmap-dir", default="/gpfs/commons/groups/sanjana_lab/Cas13/"
                    "MASLD_library_design/data/external/roadmap_liver_E066")
    ap.add_argument("--abc-bed", default="/gpfs/commons/groups/sanjana_lab/Cas13/"
                    "MASLD_library_design/data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz",
                    help="hg19 ABC liver enhancer BED (universe-independent ABC which-gene overlap)")
    ap.add_argument("--src-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--out", default=os.path.join(SF, "ag_mechanism_profile.tsv"))
    ap.add_argument("--key-file", default=None)
    ap.add_argument("--seq-length", default="SEQUENCE_LENGTH_1MB")
    ap.add_argument("--limit", type=int, default=0, help="first N regulatory leads (smoke)")
    ap.add_argument("--anchors-only", action="store_true")
    ap.add_argument("--skip-coding", action="store_true")
    ap.add_argument("--sleep", type=float, default=0.2)
    ap.add_argument("--max-retries", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true", help="build inputs, no API call")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ---- key ----
    api_key, key_src = load_api_key(args.key_file)
    if api_key is None and not args.dry_run:
        log("NO API KEY. export ALPHA_GENOME_API_KEY=<key> or write ~/.alphagenome_key "
            "(chmod 600). Key never logged.")
        sys.exit(3)
    if api_key is not None:
        log(f"API key loaded from {key_src} (len={len(api_key)}; not shown)")

    # ---- coloc status ----
    coloc = {}
    if os.path.exists(args.substrate):
        for r in csv.DictReader(open(args.substrate), delimiter="\t"):
            v = str(r.get("variant_id_hg19", "")).replace("chr", "")
            coloc[v] = ("eQTL-present" if str(r.get("colocalizes", "")).upper() == "TRUE"
                        else "eQTL-absent")

    # ---- regulatory leads ----
    leads = []
    if not args.anchors_only:
        for r in csv.DictReader(open(args.leads), delimiter="\t"):
            leads.append({
                "variant_id": r["variant_id_hg19"].replace("chr", ""),
                "variant_id_hg19": r["variant_id_hg19"].replace("chr", ""),
                "locus_id": r.get("locus_id", ""),
                "chr": str(r["chr"]).replace("chr", ""),
                "pos_hg38": int(r["pos_hg38"]),
                "hg38_ref": r.get("ref_hg38", ""), "hg38_alt": r.get("alt_hg38", ""),
                "gene": r.get("gene_substrate", ""), "ensembl": r.get("ensembl_substrate", ""),
                "lead_pip": r.get("lead_pip", ""), "var_class": "regulatory",
                "anchor_name": "",
            })
        if args.limit > 0:
            leads = leads[:args.limit]
    log(f"regulatory leads: {len(leads)}")

    # ---- Borzoi best-gene (for ag_contact_agrees_borzoi) ----
    borzoi_best = {}
    if os.path.exists(args.nominations):
        for r in csv.DictReader(open(args.nominations), delimiter="\t"):
            v = str(r.get("variant_id_hg19", "")).replace("chr", "")
            borzoi_best[v] = (str(r.get("best_gene", "")), strip_ver(str(r.get("best_gene_ensembl", ""))))

    # ---- coding effectors (42 rows; 30 in-target eQTL-absent) ----
    coding = []
    if not args.anchors_only and not args.skip_coding and os.path.exists(args.coding):
        for r in csv.DictReader(open(args.coding), delimiter="\t"):
            vhg38 = str(r.get("variant_hg38", "")).replace("chr", "")
            parts = vhg38.split(":")
            if len(parts) < 4:
                continue
            coding.append({
                "variant_id": str(r.get("variant", "")).replace("chr", ""),
                "variant_id_hg19": str(r.get("variant", "")).replace("chr", ""),
                "locus_id": vhg38,
                "chr": parts[0], "pos_hg38": int(parts[1]),
                "hg38_ref": parts[2], "hg38_alt": parts[3],
                "gene": r.get("gene", ""), "ensembl": "",
                "protein_variant": r.get("protein_variant", ""),
                "lead_pip": "", "var_class": "coding", "anchor_name": "",
            })
    log(f"coding effectors: {len(coding)}")

    # ---- anchors ----
    lead_vids = {l["variant_id"] for l in leads}
    anchors = [dict(a) for a in ANCHORS if a["variant_id"] not in lead_vids]
    for a in anchors:
        a["var_class"] = "anchor"
    work = leads + coding + anchors
    log(f"total variants to score: {len(work)} "
        f"({len(leads)} regulatory + {len(coding)} coding + {len(anchors)} anchors)")

    # ---- candidate TSS index + GTF resolution (reuse 64b) ----
    cand_idx = build_candidate_index(args.candidates)
    # inject anchor competitor TSS so the contact which-gene test is discriminating
    for a in anchors:
        comp = ANCHOR_CANDIDATES.get(a.get("anchor_name", ""), [])
        cand_idx[a["variant_id"]] = [
            {"gene": g, "ensembl": e, "biotype": "protein_coding",
             "tss_dist": "", "is_best": "", "magnitude_percentile": ""}
            for (g, e) in comp]
    b64 = load_64b(args.src_dir)
    want_ens, want_sym = set(), set()
    for v, cands in cand_idx.items():
        for c in cands:
            if c["ensembl"]:
                want_ens.add(strip_ver(c["ensembl"]))
            if c["gene"]:
                want_sym.add(c["gene"])
    for w in work:
        if w.get("ensembl"):
            want_ens.add(strip_ver(w["ensembl"]))
        if w.get("gene"):
            want_sym.add(w["gene"])
    log(f"resolving TSS for {len(want_ens)} ensembl / {len(want_sym)} symbols from GTF ...")
    genes, sym2sid = b64.parse_gtf_genes(args.gtf, want_ens, want_sym)

    def resolve_tss(ensembl, symbol, chrom):
        rec, how = b64.resolve_gene(genes, sym2sid, ensembl, symbol,
                                    chrom if chrom.startswith("chr") else "chr" + chrom)
        if rec is None or rec.get("gene_start") is None:
            return None, None, how
        tss1, strand = b64.gene_tss(rec)
        return tss1, strand, how

    # ---- corroborators ----
    abc = load_abc()
    abc_bed = load_abc_bed(args.abc_bed)
    ccre = load_ccre()
    cbp = load_chrombpnet()
    mbr = load_motifbreakr()
    roadmap = load_roadmap_beds(args.roadmap_dir)
    log(f"corroborators: ABC precomputed {len(abc)} vars + ABC-bed {sum(len(v) for v in abc_bed.values())} "
        f"enhancers/{len(abc_bed)} chr, cCRE {len(ccre)}, ChromBPNet {len(cbp)}, "
        f"motifbreakR {len(mbr)}, roadmap heads {sorted(roadmap.keys())}")

    if args.dry_run:
        n_ok = sum(1 for w in work if w.get("hg38_ref") and w.get("hg38_alt"))
        log(f"DRY-RUN: {n_ok}/{len(work)} variants carry hg38 alleles; "
            f"{sum(len(cand_idx.get(w['variant_id'], [])) for w in leads)} candidate TSS across leads. "
            f"No API call.")
        sys.exit(0)

    # ---- model ----
    from alphagenome.data import genome
    from alphagenome.models import dna_client, dna_output, variant_scorers
    import grpc
    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[args.seq_length]
    R = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    scalar_scorers = [R[k] for k in SCALAR_SCORER_KEYS if k in R]
    model = dna_client.create(api_key)
    log(f"dna_client created; seq_len={seq_len:,}; {len(scalar_scorers)} scalar scorers + CONTACT_MAPS")

    def _retry(fn, tag):
        last = None
        for attempt in range(args.max_retries):
            try:
                return fn()
            except grpc.RpcError as e:  # noqa
                last = e
                code = e.code()
                if code not in (grpc.StatusCode.UNAVAILABLE, grpc.StatusCode.RESOURCE_EXHAUSTED,
                                grpc.StatusCode.DEADLINE_EXCEEDED, grpc.StatusCode.INTERNAL):
                    raise
                wait = min(45.0, 2.0 ** attempt)
                log(f"  transient {code} on {tag} (attempt {attempt+1}); backoff {wait:.0f}s")
                time.sleep(wait)
        raise last

    def score_fingerprint(chrom, pos1, ref, alt, vid):
        v = genome.Variant(chromosome=chrom, position=pos1, reference_bases=ref,
                           alternate_bases=alt, name=vid)
        iv = v.reference_interval.resize(seq_len)
        sc = _retry(lambda: model.score_variant(
            interval=iv, variant=v, variant_scorers=scalar_scorers,
            organism=dna_client.Organism.HOMO_SAPIENS), f"fingerprint {vid}")
        return variant_scorers.tidy_scores(sc)

    def predict_contact(chrom, pos1, ref, alt, vid):
        v = genome.Variant(chromosome=chrom, position=pos1, reference_bases=ref,
                           alternate_bases=alt, name=vid)
        iv = v.reference_interval.resize(seq_len)
        vout = _retry(lambda: model.predict_variant(
            interval=iv, variant=v, requested_outputs=[dna_output.OutputType.CONTACT_MAPS],
            ontology_terms=None), f"contact {vid}")
        return vout, iv

    # HepG2 contact-track index (resolved lazily from the first contact metadata)
    hepg2_idx = None

    rows = []
    n_fp = n_ct = n_err = 0
    dbg = True
    for k, w in enumerate(work):
        vid = w["variant_id"]
        chrom_raw = str(w["chr"]).replace("chr", "")
        chrom = "chr" + chrom_raw
        pos1 = int(w["pos_hg38"])
        ref, alt = str(w.get("hg38_ref", "")).upper(), str(w.get("hg38_alt", "")).upper()
        var_class = w.get("var_class", "regulatory")
        is_anchor = var_class == "anchor"
        target_gene = w.get("gene", "")
        target_ens = w.get("ensembl", "")
        note = ""
        if not ref or not alt:
            note = "missing_hg38_allele"

        row = {"variant_id": vid, "variant_id_hg19": w.get("variant_id_hg19", vid),
               "var_class": var_class, "chr": chrom_raw, "pos_hg38": pos1,
               "hg38_ref": ref, "hg38_alt": alt,
               "target_gene": target_gene, "target_ensembl": target_ens,
               "coloc_status": coloc.get(vid, "unknown"),
               "in_target": str(coloc.get(vid, "") == "eQTL-absent").upper(),
               "lead_pip": w.get("lead_pip", ""),
               "protein_variant": w.get("protein_variant", ""),
               "is_anchor": str(is_anchor).upper(), "anchor_name": w.get("anchor_name", ""),
               "note": note}

        if not ref or not alt:
            rows.append(row)
            continue

        # ---- C1 fingerprint ----
        try:
            df = score_fingerprint(chrom, pos1, ref, alt, vid)
            if dbg:
                log(f"  fingerprint tidy cols OK ({len(df)} rows)")
                dbg = False
            fp = extract_fingerprint(df, target_ens or None, target_gene or None)
            row.update({k2: v2 for k2, v2 in fp.items()})
            n_fp += 1
        except Exception as e:  # noqa
            row["note"] = (note + ";" if note else "") + f"fingerprint_err:{type(e).__name__}:{e}"
            rows.append(row)
            n_err += 1
            time.sleep(args.sleep)
            continue

        # ---- C2 contact which-gene (regulatory leads + anchors; coding = self) ----
        cands = cand_idx.get(vid, [])
        if var_class == "coding":
            row.update(ag_contact_gene=target_gene, ag_contact_gene_ensembl="",
                       ag_contact_gene_tss_dist="", contact_rewire_abs="",
                       contact_track_used="n/a(coding)", n_candidate_tss=0,
                       nearest_tss_gene=target_gene, ag_contact_overturns_nearest="FALSE",
                       ag_contact_agrees_hepg2="")
        elif len(cands) == 0:
            row.update(ag_contact_gene="", ag_contact_gene_ensembl="",
                       ag_contact_gene_tss_dist="", contact_rewire_abs="",
                       contact_track_used="", n_candidate_tss=0, nearest_tss_gene="",
                       ag_contact_overturns_nearest="", ag_contact_agrees_hepg2="")
            row["note"] = (row.get("note", "") + ";" if row.get("note") else "") + "no_candidates"
        else:
            try:
                vout, iv = predict_contact(chrom, pos1, ref, alt, vid)
                ref_cm = vout.reference.contact_maps
                alt_cm = vout.alternate.contact_maps
                cvals = np.asarray(ref_cm.values, dtype=np.float32)   # (nb,nb,ntracks)
                avals = np.asarray(alt_cm.values, dtype=np.float32)
                nb = cvals.shape[0]
                binsize = iv.width / nb
                cstart = ref_cm.interval.start if hasattr(ref_cm, "interval") else iv.start
                if hepg2_idx is None and hasattr(ref_cm, "metadata"):
                    mdc = ref_cm.metadata
                    hm = liver_mask(mdc)
                    hepg2_idx = int(np.where(hm)[0][0]) if hm.sum() > 0 else -1
                    log(f"  contact tracks={cvals.shape[2]} binsize={binsize:.0f}bp "
                        f"HepG2_track_idx={hepg2_idx}")
                ref_avg = np.nanmean(cvals, axis=2)     # tissue-average (primary)
                alt_avg = np.nanmean(avals, axis=2)
                ref_hep = cvals[:, :, hepg2_idx] if (hepg2_idx is not None and hepg2_idx >= 0) else None
                vb = int((pos1 - cstart) // binsize)
                vb = min(max(vb, 0), nb - 1)

                best = None       # (contact, cand, tss_bin, tss1)
                best_hep = None
                nearest = None    # (abs_tss_dist, cand)
                for c in cands:
                    tss1, strand, how = resolve_tss(c["ensembl"], c["gene"], chrom_raw)
                    if tss1 is None:
                        continue
                    tb = int((tss1 - 1 - cstart) // binsize)
                    if tb < 0 or tb >= nb:
                        continue
                    contact = float(ref_avg[vb, tb])
                    try:
                        td = abs(int(float(c.get("tss_dist", "nan"))))
                    except Exception:
                        td = abs(pos1 - tss1)
                    if best is None or contact > best[0]:
                        best = (contact, c, tb, tss1)
                    if ref_hep is not None:
                        ch = float(ref_hep[vb, tb])
                        if best_hep is None or ch > best_hep[0]:
                            best_hep = (ch, c)
                    if nearest is None or td < nearest[0]:
                        nearest = (td, c)

                if best is None:
                    row.update(ag_contact_gene="", ag_contact_gene_ensembl="",
                               ag_contact_gene_tss_dist="", contact_rewire_abs="",
                               contact_track_used="tissue_avg(28)", n_candidate_tss=len(cands),
                               nearest_tss_gene=(nearest[1]["gene"] if nearest else ""),
                               ag_contact_overturns_nearest="", ag_contact_agrees_hepg2="")
                    row["note"] = (row.get("note", "") + ";" if row.get("note") else "") + "no_candidate_in_binspan"
                else:
                    contact, cbest, tb, tss1 = best
                    rewire = float(abs(alt_avg[vb, tb] - ref_avg[vb, tb]))
                    nearest_gene = nearest[1]["gene"] if nearest else ""
                    overturns = (cbest["gene"] != nearest_gene) if nearest_gene else False
                    hep_gene = best_hep[1]["gene"] if best_hep else ""
                    row.update(
                        ag_contact_gene=cbest["gene"],
                        ag_contact_gene_ensembl=strip_ver(cbest["ensembl"]),
                        ag_contact_gene_tss_dist=cbest.get("tss_dist", ""),
                        contact_rewire_abs=f"{rewire:.6f}",
                        contact_track_used="tissue_avg(28)",
                        n_candidate_tss=len(cands),
                        nearest_tss_gene=nearest_gene,
                        ag_contact_overturns_nearest=str(overturns).upper(),
                        ag_contact_agrees_hepg2=str(bool(hep_gene) and hep_gene == cbest["gene"]).upper(),
                        contact_diag_oe=f"{float(ref_avg[vb, vb]):.4f}",   # O/E self-bin (neg => distance-normalised map)
                    )
                    # agreement vs Borzoi best gene
                    bz_gene, bz_ens = borzoi_best.get(vid, ("", ""))
                    agrees_bz = bool(cbest["gene"]) and (
                        (bz_ens and strip_ver(cbest["ensembl"]) == bz_ens) or
                        (bz_gene and cbest["gene"] == bz_gene))
                    row["ag_contact_agrees_borzoi"] = str(agrees_bz).upper()
                    row["borzoi_best_gene"] = bz_gene
                n_ct += 1
            except Exception as e:  # noqa
                row["note"] = (row.get("note", "") + ";" if row.get("note") else "") + \
                    f"contact_err:{type(e).__name__}:{e}"

        rows.append(row)
        if is_anchor or k % 15 == 0:
            log(f"[{k+1}/{len(work)}] {vid} {target_gene} {var_class} "
                f"contact_gene={row.get('ag_contact_gene','')} "
                f"q_rna={row.get('q_rna_gene','')} q_tf={row.get('q_chip_tf','')} "
                f"splice={max(_absnan(row.get('q_splice_junctions')), _absnan(row.get('q_splice_sites')), _absnan(row.get('q_splice_site_usage'))):.3f} "
                f"{'['+row['note']+']' if row.get('note') else ''}")
        time.sleep(args.sleep)

    df_out = pd.DataFrame(rows)

    # ---- C2 cohort contact percentile (regulatory only) ----
    cr = pd.to_numeric(df_out.get("contact_rewire_abs", pd.Series(dtype=float)), errors="coerce")
    reg_mask = df_out["var_class"].isin(["regulatory", "anchor"]) & cr.notna()
    df_out["contact_pctile"] = ""
    if reg_mask.sum() > 1:
        ranks = cr[reg_mask].rank(pct=True) * 100.0
        df_out.loc[reg_mask, "contact_pctile"] = ranks.round(1).astype(str)

    # ---- C3 mechanism class + C4 corroboration (row-wise post-pass) ----
    mech_cols = {"mechanism_class": [], "mechanism_class_secondary": [], "mechanism_confidence": [],
                 "splice_score": [], "promoter_score": [], "enhancer_score": [], "tf_score": [],
                 "ag_element_class": [], "roadmap_state": [], "roadmap_h3k27ac": [], "roadmap_h3k4me3": [],
                 "ccre_class": [], "ccre_gene": [], "ag_contact_agrees_abc": [], "abc_target_gene": [],
                 "ag_tf_chip_disrupted": [], "chrombpnet_disrupted_tf": [], "motifbreakr_disrupted_tf": []}
    for _, r in df_out.iterrows():
        vid = r["variant_id"]
        fp = {c: r.get(c) for c in ("q_splice_junctions", "q_splice_sites", "q_splice_site_usage",
                                    "q_cage", "q_procap", "q_h3k4me3", "q_atac", "q_dnase",
                                    "q_h3k27ac", "q_h3k4me1", "q_chip_tf")}
        try:
            tss_dist = abs(int(float(r.get("ag_contact_gene_tss_dist", "") or "nan")))
        except Exception:
            tss_dist = None
        try:
            cp = float(r.get("contact_pctile", "") or "nan")
        except Exception:
            cp = float("nan")
        contact_high = (cp == cp) and cp >= CONTACT_HIGH_PCTILE
        # ---- C4 corroborators FIRST (feed the class rule) ----
        chromhmm_state = roadmap_hit(roadmap, "chromHMM", str(r["chr"]), int(r["pos_hg38"]))
        cbp_tf = cbp.get(vid, "")
        mbr_tfs = mbr.get(vid, set())
        top_tf = str(r.get("top_tf", ""))
        tf_disrupted = bool(cbp_tf) or bool(mbr_tfs)
        tf_corroborated = tf_disrupted and bool(top_tf)
        # ---- C3 mechanism class (corroboration-aware) ----
        if r["var_class"] == "coding":
            prim, sec, conf = "protein-altering", "", ""
            ss = ps = es = ts = ""
        else:
            prim, sec, conf, sc = assign_mechanism_class(
                fp, tss_dist, contact_high, tf_corroborated=tf_corroborated,
                chromhmm_state=chromhmm_state)
            ss, ps, es, ts = (round(sc["splice-altering"], 4), round(sc["promoter/TSS"], 4),
                              round(sc["enhancer-disrupting"], 4), round(sc["TF-footprint-breaking"], 4))
        mech_cols["mechanism_class"].append(prim)
        mech_cols["mechanism_class_secondary"].append(sec)
        mech_cols["mechanism_confidence"].append(conf)
        mech_cols["splice_score"].append(ss)
        mech_cols["promoter_score"].append(ps)
        mech_cols["enhancer_score"].append(es)
        mech_cols["tf_score"].append(ts)
        elt = {"promoter/TSS": "promoter", "enhancer-disrupting": "enhancer",
               "TF-footprint-breaking": "TF_bound", "splice-altering": "splice",
               "3D-contact-rewiring": "loop_anchor"}.get(prim, "") if r["var_class"] != "coding" else ""
        mech_cols["ag_element_class"].append(elt)
        # Roadmap E066 + SCREEN cCRE (measured, independent)
        _st = str(chromhmm_state).strip()
        mech_cols["roadmap_state"].append(f"{_st}:{CHROMHMM_NAMES[_st]}" if _st in CHROMHMM_NAMES else _st)
        mech_cols["roadmap_h3k27ac"].append("TRUE" if roadmap_hit(roadmap, "H3K27ac", str(r["chr"]), int(r["pos_hg38"])) else "FALSE")
        mech_cols["roadmap_h3k4me3"].append("TRUE" if roadmap_hit(roadmap, "H3K4me3", str(r["chr"]), int(r["pos_hg38"])) else "FALSE")
        cc = ccre.get(vid, ("", ""))
        mech_cols["ccre_class"].append(cc[0])
        mech_cols["ccre_gene"].append(cc[1])
        # ABC which-gene corroboration (independent liver AvgHiC): precomputed exact-variant
        # join UNION direct hg19 enhancer-BED overlap (universe-independent, wider coverage).
        abc_genes = set(abc.get(vid, set()))
        try:
            _pos_hg19 = int(str(r.get("variant_id_hg19", "")).split(":")[1])
            abc_genes |= abc_bed_genes(abc_bed, str(r["chr"]), _pos_hg19)
        except (ValueError, IndexError):
            pass
        cg = str(r.get("ag_contact_gene", ""))
        mech_cols["abc_target_gene"].append(";".join(sorted(abc_genes)) if abc_genes else "")
        mech_cols["ag_contact_agrees_abc"].append(
            "TRUE" if (cg and cg in abc_genes) else ("FALSE" if abc_genes else ""))
        mech_cols["ag_tf_chip_disrupted"].append(str(tf_corroborated).upper())
        mech_cols["chrombpnet_disrupted_tf"].append(cbp_tf)
        mech_cols["motifbreakr_disrupted_tf"].append(";".join(sorted(mbr_tfs)) if mbr_tfs else "")
    for c, vals in mech_cols.items():
        df_out[c] = vals

    # ---- column order + write ----
    front = ["variant_id", "variant_id_hg19", "var_class", "coloc_status", "in_target",
             "chr", "pos_hg38", "hg38_ref", "hg38_alt", "target_gene", "target_ensembl",
             "protein_variant", "lead_pip", "is_anchor", "anchor_name",
             "mechanism_class", "mechanism_class_secondary", "mechanism_confidence",
             "ag_contact_gene", "ag_contact_gene_ensembl", "ag_contact_gene_tss_dist",
             "nearest_tss_gene", "ag_contact_overturns_nearest", "borzoi_best_gene",
             "ag_contact_agrees_borzoi", "ag_contact_agrees_abc", "ag_contact_agrees_hepg2",
             "abc_target_gene", "contact_rewire_abs", "contact_diag_oe", "contact_pctile",
             "contact_track_used", "n_candidate_tss", "ag_element_class",
             "q_rna_gene", "q_polya", "q_splice_junctions", "q_splice_sites", "q_splice_site_usage",
             "q_atac", "q_dnase", "q_chip_histone", "q_chip_tf", "q_cage", "q_procap",
             "q_h3k4me3", "q_h3k27ac", "q_h3k4me1",
             "splice_score", "promoter_score", "enhancer_score", "tf_score",
             "top_tf", "top_tf_quantile", "top_histone_mark", "top_histone_quantile",
             "ag_tf_chip_disrupted", "chrombpnet_disrupted_tf", "motifbreakr_disrupted_tf",
             "roadmap_state", "roadmap_h3k27ac", "roadmap_h3k4me3", "ccre_class", "ccre_gene",
             "n_liver_rna", "note"]
    cols = [c for c in front if c in df_out.columns] + [c for c in df_out.columns if c not in front]
    df_out = df_out[cols]
    df_out.to_csv(args.out, sep="\t", index=False)
    log(f"WROTE {args.out} ({len(df_out)} rows, {len(cols)} cols)")

    # ---- README ----
    write_readme(df_out, args.out)

    # ---- anchor verdict ----
    log("================= ANCHOR CHECK =================")
    a_sort1 = df_out[df_out["anchor_name"] == "SORT1_posctrl"]
    if len(a_sort1):
        r = a_sort1.iloc[0]
        ok = str(r.get("ag_contact_gene", "")) == "SORT1"
        log(f"SORT1: ag_contact_gene={r.get('ag_contact_gene','')} "
            f"(expect SORT1) -> {'PASS' if ok else 'CHECK'}; "
            f"top_tf={r.get('top_tf','')} q={r.get('top_tf_quantile','')}")
    a_hsd = df_out[df_out["anchor_name"] == "HSD17B13_splice"]
    if len(a_hsd):
        r = a_hsd.iloc[0]
        sq = max(_absnan(r.get("q_splice_junctions")), _absnan(r.get("q_splice_sites")),
                 _absnan(r.get("q_splice_site_usage")))
        log(f"HSD17B13: splice |quantile|={sq:.4f} (expect ~1.0) -> "
            f"{'PASS' if sq >= 0.9 else 'CHECK'}; mechanism={r.get('mechanism_class','')}")
    log("===============================================")
    log(f"fingerprints scored: {n_fp}; contact predicted: {n_ct}; errors: {n_err}")


def write_readme(df, out_tsv):
    md_path = out_tsv.replace(".tsv", ".README.md")
    reg = df[df["var_class"] == "regulatory"]
    cod = df[df["var_class"] == "coding"]
    mc = reg["mechanism_class"].value_counts().to_dict() if "mechanism_class" in reg else {}
    with_contact = reg[reg.get("ag_contact_gene", "").astype(str).str.len() > 0]
    n_overturn = int((reg.get("ag_contact_overturns_nearest", pd.Series(dtype=str)) == "TRUE").sum())
    bz_ov = reg[reg.get("ag_contact_agrees_borzoi", "").isin(["TRUE", "FALSE"])]
    bz_rate = f"{int((bz_ov['ag_contact_agrees_borzoi']=='TRUE').sum())}/{len(bz_ov)}" if len(bz_ov) else "0/0"
    abc_ov = reg[reg.get("ag_contact_agrees_abc", "").isin(["TRUE", "FALSE"])]
    abc_rate = (f"{int((abc_ov['ag_contact_agrees_abc']=='TRUE').sum())}/{len(abc_ov)}"
                if len(abc_ov) else "0/0 (no ABC overlap)")
    md = f"""# AlphaGenome multimodal mechanism profile (Phase-4, script 72)

**Mechanism-CLASS / magnitude ONLY. NOT a direction claim, NOT a scored convergence
channel.** The eQTL-direction gate FAILED for both sequence models (Borzoi auROC 0.559 /
AlphaGenome 0.561, ns), so signed quantiles are recorded but never asserted as up/down.
Folded into src/70 as Axis-5 ANNOTATIONS; absent from convergence_evidence.csv inputs.

## Rows
- {len(reg)} regulatory eQTL-absent leads (fingerprint + mechanism-class + contact loop-partner)
- {len(cod)} coding effectors (fingerprint + protein-altering; {int((cod['in_target']=='TRUE').sum())} in-target eQTL-absent)
- {int((df['var_class']=='anchor').sum())} anchors (SORT1 / HSD17B13 / PNPLA3)

## Components
- **C1 fingerprint (SOLID)**: one score_variant/variant, 11 liver-masked scalar heads ->
  **MEAN signed quantile** over liver tracks (calibrated + multiplicity-robust; the earlier
  signed-EXTREME saturated at ~1.0 for any variant when a head has 100s of liver tracks).
  `splice_junctions` is EXCLUDED from the class rule (degenerate null -> saturates regardless).
- **C3 mechanism class (SOLID, corroboration-aware)**: argmax over discriminative MEAN |quantiles|
  {{splice_sites/usage, promoter marks, enhancer marks, TF}}; TF-footprint asserted ONLY when an
  independent footprint caller (ChromBPNet / motifbreakR) corroborates; MEASURED Roadmap chromHMM
  state disambiguates promoter vs enhancer. Distribution: {mc}
- **C4 corroboration (SOLID)**: measured, independent -> agreement FLAGS only (Roadmap E066 chromHMM +
  H3K27ac/H3K4me3, SCREEN cCRE, ChromBPNet, motifbreakR).
- **C2 contact which-gene (NEGATIVE / exploratory only -- see below).**

## C2 contact which-gene: HONEST NEGATIVE RESULT (do NOT use to supersede nominations)
AlphaGenome's contact map is **O/E-normalised** (the self-bin O/E is NEGATIVE: SORT1 -0.55,
KLHL8 -0.55, FCGRT -0.30, ...), i.e. distance-corrected. Argmaxing candidate-TSS contact from
the variant bin is therefore a **long-range LOOP-PARTNER detector, not an effector-gene caller**:
it systematically over-calls DISTAL genes ({n_overturn}/{len(with_contact)} "overturn" nearest-TSS,
mostly toward distal / non-coding TSS) while the proximal effector sits at ~0 or negative O/E.
- AG-contact vs Borzoi-magnitude which-gene: **{bz_rate}** (~chance).
- AG-contact vs ABC (MEASURED liver AvgHiC): **{abc_rate}** -- on the one lead where ABC has an
  opinion (KLHL8) ABC calls KLHL8 (the nearest gene) and AG-contact calls a 235 kb-distal gene.
- **PASS on SORT1**: contact recovers the true distal effector SORT1 (123 kb) over nearest PSRC1
  (9 kb) / CELSR2 -- a genuine strong enhancer-promoter loop. So contact WORKS for real long-range
  loops but is NOT a reliable blanket which-gene method; it is retained as an EXPLORATORY
  "candidate long-range loop partner" annotation, class-only, and does **NOT** replace the
  nearest-TSS / Borzoi which-gene call in src/70. A validated version would need ABC-style
  activity x contact weighting (out of scope; the raw O/E argmax is reported as-is, honestly).

## Guardrails
- Corroborators are agreement flags, never summed. AG-contact scored vs Borzoi + ABC (independent),
  never vs AG's own expression head (circularity discipline).
- AG outputs must not train other models (license); research-only. Needs orthogonal wet-lab / MPRA.

_Generated by src/72_ag_multimodal_mechanism.py (additive; not wired into 46d/78/27a)._
"""
    with open(md_path, "w") as fh:
        fh.write(md)
    log(f"WROTE {md_path}")


if __name__ == "__main__":
    main()
