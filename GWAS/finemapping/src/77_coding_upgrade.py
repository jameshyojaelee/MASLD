#!/usr/bin/env python3
"""
77_coding_upgrade.py  --  Axis-1 CODING-ARM UPGRADE (Phase-5 Track 3).

STRICTLY ADDITIVE. Does NOT touch 46d/78/27a or the 60-73 seqfunc scripts.
Reads the finished coding effectors (coding_hardening.tsv, 42 rows) + the
coding_rescue.tsv superset, and adds two independent channels the current
LoF/pathogenicity-oriented consensus (AlphaMissense + ESM1b + popEVE) misses:

  1. GAIN- vs LOSS-of-function DIRECTION  --  LoGoFunc (Stein et al. 2023,
     Genome Medicine s13073-023-01261-9). Precomputed genome-wide (hg38) GoF /
     LoF / neutral probabilities for 71.3M canonical missense variants
     (Zenodo 10126185, tabix-indexed). We query the ~42 effector coordinates
     and take argmax(neutral, GOF, LOF) as the raw LoGoFunc call.

     WHY LoGoFunc and not PreMode: PreMode (Shen lab, Nat Commun 2025) is the
     purpose-built mode-of-action model, BUT it ships precomputed scores for
     only 10 genes (BRAF/RET/TP53/KCNJ11/CACNA1A/SCN5A/SCN2A/ABCC8/PTEN);
     PNPLA3 and our other effectors are NOT precomputed. Running PreMode on a
     new gene needs per-gene AlphaFold2 structure + MSA + GPU inference set up
     from its HuggingFace weights -- for 42 genes that is well beyond a
     reasonable additive effort. LoGoFunc is genome-wide, turnkey, and gives a
     directional GoF/LoF call for every canonical missense. (Documented so a
     future run can swap in PreMode for a targeted gene if warranted.)

  2. EVE  --  Frazer et al. 2021, Nature s41586-021-04043-8. An evolution-ONLY
     (deep generative model of an MSA) severity + Benign/Uncertain/Pathogenic
     class, precomputed for ~3,219 disease genes (evemodel.org, GRCh38 per-
     protein VCFs). This is an INDEPENDENT predictor from the existing
     consensus: AlphaMissense uses supervised + structure, ESM1b is a protein
     LLM, popEVE is EVE+ESM re-calibrated to human population. EVE-alone adds a
     pure-evolution vote and covers effectors AM missed (bonus: PARVB).

MERE NOMINATION / MECHANISM-CLASS ONLY. These columns are protein-model
severity + direction for the eQTL-blind coding effectors; they are NEVER wired
into the convergence atlas or a scored channel (would be circular vs COLOC /
epigenomic).

OUTPUT
  coding_hardening_v2.tsv  = coding_hardening.tsv + new columns:
     logofunc_neutral, logofunc_gof, logofunc_lof, logofunc_call, logofunc_note
     gof_lof_call    -- FINAL direction: curated_literature anchor if present
                        (experimentally established), else LoGoFunc argmax
     gof_lof_source  -- curated_literature / LoGoFunc / none
     eve_score, eve_class, eve_note
  coding_hardening_v2.README.md

RUN (rnaseq env; pysam/tabix present, htslib on PATH)
  micromamba run -n rnaseq python GWAS/finemapping/src/77_coding_upgrade.py
Optional: --eve-dir <dir with per-protein EVE VCFs> (default data/external/eve/vcf).
Re-run after EVE lands to fill eve_* columns (LoGoFunc columns are stable).
"""
import os
import re
import csv
import glob
import argparse
import subprocess

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SEQFUNC = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
CODING_HARDENING = os.path.join(SEQFUNC, "coding_hardening.tsv")
CODING_RESCUE = os.path.join(SEQFUNC, "coding_rescue.tsv")
LOGOFUNC_GZ = os.path.join(
    ROOT, "data/external/logofunc/LoGoFuncVotingEnsemble_preds_final.csv.gz"
)
OUT_TSV = os.path.join(SEQFUNC, "coding_hardening_v2.tsv")
OUT_README = os.path.join(SEQFUNC, "coding_hardening_v2.README.md")

# htslib tabix binary (module load samtools / htslib provides it)
TABIX = "tabix"
for cand in (
    "/nfs/sw/easybuild/software/htslib/1.23.1/bin/tabix",
):
    if os.path.exists(cand):
        TABIX = cand
        break


def parse_variant_hg38(v):
    """'22:43928847:C:G' -> ('22', 43928847, 'C', 'G')."""
    parts = v.split(":")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    chrom = chrom.replace("chr", "")
    try:
        pos = int(pos)
    except ValueError:
        return None
    return chrom, pos, ref, alt


# ---------------------------------------------------------------------------
# LoGoFunc: tabix query by hg38 coordinate, match REF+ALT, argmax the 3 probs.
# File columns: #CHROM POS REF ALT ID prediction LoGoFunc_neutral GOF LOF
# ---------------------------------------------------------------------------
def logofunc_lookup(chrom, pos, ref, alt):
    try:
        out = subprocess.run(
            [TABIX, LOGOFUNC_GZ, f"{chrom}:{pos}-{pos}"],
            capture_output=True, text=True, timeout=60,
        ).stdout
    except Exception as e:
        return dict(note=f"tabix_error:{e}")
    hit = None
    for line in out.strip().splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        c, p, r, a = f[0], f[1], f[2], f[3]
        if p == str(pos) and r == ref and a == alt:
            hit = f
            break
    if hit is None:
        # position present but this REF/ALT not among LoGoFunc's canonical
        # missense alleles (e.g. non-canonical transcript / not missense in
        # canonical) -> absent
        return dict(note="absent_from_logofunc")
    neutral = float(hit[6])
    gof = float(hit[7])
    lof = float(hit[8])
    call_map = {"neutral": "Neutral", "gof": "GOF", "lof": "LOF"}
    argmax = max((neutral, "neutral"), (gof, "gof"), (lof, "lof"))[1]
    return dict(
        neutral=neutral, gof=gof, lof=lof,
        call=call_map[argmax],
        raw_prediction=hit[5],
        note="covered",
    )


# ---------------------------------------------------------------------------
# EVE: per-protein GRCh38 VCFs (evemodel.org bulk) are UniProt-entry-named
# (PLPL3_HUMAN.vcf, not PNPLA3_HUMAN.vcf), so filename->gene matching is
# unreliable. We instead merge all per-protein VCFs into one bgzip+tabix file
# (see 77_build_eve_merged.sh) and look variants up by GRCh38 coordinate --
# gene-agnostic, exactly the Ensembl-VEP EVE merge pattern. INFO carries
# EVE_SCORE + EVE_CLASS (Benign/Uncertain/Pathogenic, Class75 default).
# Optional -- skipped gracefully if the merged EVE file is not on disk yet.
# ---------------------------------------------------------------------------
EVE_MERGED = os.path.join(ROOT, "data/external/eve/eve_merged.vcf.gz")


def _parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def eve_lookup(chrom, pos, ref, alt, protein_variant):
    try:
        out = subprocess.run(
            [TABIX, EVE_MERGED, f"{chrom}:{pos}-{pos}"],
            capture_output=True, text=True, timeout=60,
        ).stdout
    except Exception as e:
        return dict(note=f"eve_tabix_error:{e}")
    coord_hit = None
    aa_hit = None
    for line in out.strip().splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        vc, vp, _id, vr, va = f[0].replace("chr", ""), f[1], f[2], f[3], f[4]
        d = _parse_info(f[7])
        if vc == str(chrom) and vp == str(pos) and vr == ref and va == alt:
            coord_hit = d
            break
        aa = (d.get("EVE_PROT_MUT") or d.get("MUT") or d.get("wt_aa", ""))
        if protein_variant and aa and protein_variant in aa:
            aa_hit = d
    d = coord_hit or aa_hit
    if d is None:
        note = ("variant_not_in_eve" if out.strip()
                else "gene_not_in_eve_or_uncovered_pos")
        return dict(note=note)
    return dict(score=d.get("EVE_SCORE", ""), cls=d.get("EVE_CLASS", ""),
                note=("covered_coord" if coord_hit else "covered_aa"))


def main():
    global EVE_MERGED
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--eve-merged", default=EVE_MERGED,
        help="bgzip+tabix merged EVE GRCh38 VCF (eve_merged.vcf.gz)",
    )
    args = ap.parse_args()

    if not os.path.exists(LOGOFUNC_GZ):
        raise SystemExit(f"LoGoFunc file missing: {LOGOFUNC_GZ}")

    EVE_MERGED = args.eve_merged
    eve_available = os.path.exists(EVE_MERGED) and os.path.exists(EVE_MERGED + ".tbi")
    print(f"[eve] merged EVE file {EVE_MERGED} present -- filling eve_* columns"
          if eve_available else f"[eve] merged EVE file not found ({EVE_MERGED}) "
          "-- eve_* columns left blank (re-run when EVE lands)")

    with open(CODING_HARDENING) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    base_cols = list(rows[0].keys())
    print(f"[in] {len(rows)} coding effectors from coding_hardening.tsv")

    new_cols = [
        "logofunc_neutral", "logofunc_gof", "logofunc_lof",
        "logofunc_call", "logofunc_note",
        "gof_lof_call", "gof_lof_source",
        "eve_score", "eve_class", "eve_note",
    ]

    gof_dist = {}
    logofunc_dist = {}
    for r in rows:
        pv = parse_variant_hg38(r.get("variant_hg38", ""))
        if pv is None:
            for c in new_cols:
                r[c] = ""
            r["logofunc_note"] = "bad_coord"
            r["eve_note"] = "bad_coord"
            r["gof_lof_call"] = ""
            r["gof_lof_source"] = "none"
            continue
        chrom, pos, ref, alt = pv

        # --- LoGoFunc ---
        lf = logofunc_lookup(chrom, pos, ref, alt)
        r["logofunc_neutral"] = f"{lf['neutral']:.4f}" if "neutral" in lf else ""
        r["logofunc_gof"] = f"{lf['gof']:.4f}" if "gof" in lf else ""
        r["logofunc_lof"] = f"{lf['lof']:.4f}" if "lof" in lf else ""
        r["logofunc_call"] = lf.get("call", "")
        r["logofunc_note"] = lf.get("note", "")
        if lf.get("call"):
            logofunc_dist[lf["call"]] = logofunc_dist.get(lf["call"], 0) + 1

        # --- FINAL gof_lof_call: prefer experimentally curated anchor ---
        curated = (r.get("gof_lof_flag") or "").strip()
        if curated:
            # e.g. 'GoF_neomorph [curated_literature]' -> keep the token
            r["gof_lof_call"] = curated.replace(" [curated_literature]", "")
            r["gof_lof_source"] = "curated_literature"
        elif lf.get("call"):
            r["gof_lof_call"] = lf["call"]
            r["gof_lof_source"] = "LoGoFunc"
        else:
            r["gof_lof_call"] = ""
            r["gof_lof_source"] = "none"
        key = r["gof_lof_call"] or "unresolved"
        gof_dist[key] = gof_dist.get(key, 0) + 1

        # --- EVE ---
        if eve_available:
            ev = eve_lookup(chrom, pos, ref, alt, r.get("protein_variant", ""))
            r["eve_score"] = ev.get("score", "")
            r["eve_class"] = ev.get("cls", "")
            r["eve_note"] = ev.get("note", "")
        else:
            r["eve_score"] = ""
            r["eve_class"] = ""
            r["eve_note"] = "eve_data_not_downloaded"

    out_cols = base_cols + [c for c in new_cols if c not in base_cols]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in out_cols})
    print(f"[out] {OUT_TSV}  ({len(rows)} rows, +{len(new_cols)} cols)")

    # -- report --
    print("\n=== FINAL gof_lof_call distribution (curated-preferred) ===")
    for k, v in sorted(gof_dist.items(), key=lambda x: -x[1]):
        print(f"  {k:28s} {v}")
    print("\n=== raw LoGoFunc argmax distribution (all 42) ===")
    for k, v in sorted(logofunc_dist.items(), key=lambda x: -x[1]):
        print(f"  {k:10s} {v}")
    gof_genes = [
        f"{r['gene']}:{r['protein_variant']} "
        f"(final={r['gof_lof_call']}/{r['gof_lof_source']}, "
        f"LoGoFunc={r['logofunc_call']} "
        f"n{r['logofunc_neutral']}/g{r['logofunc_gof']}/l{r['logofunc_lof']})"
        for r in rows
        if "GoF" in (r["gof_lof_call"] or "") or r["logofunc_call"] == "GOF"
    ]
    print("\n=== GoF-flagged effectors (final GoF OR raw LoGoFunc GOF) ===")
    for g in gof_genes:
        print("  " + g)

    # PNPLA3 verbatim
    for r in rows:
        if r["gene"] == "PNPLA3" and r["protein_variant"] == "I148M":
            print("\n=== PNPLA3 I148M (verification target) ===")
            print(f"  LoGoFunc: call={r['logofunc_call']} "
                  f"neutral={r['logofunc_neutral']} GOF={r['logofunc_gof']} "
                  f"LOF={r['logofunc_lof']} ({r['logofunc_note']})")
            print(f"  final gof_lof_call={r['gof_lof_call']} "
                  f"source={r['gof_lof_source']}")
            print(f"  EVE: score={r['eve_score']} class={r['eve_class']} "
                  f"({r['eve_note']})")

    # PARVB (bonus)
    for r in rows:
        if r["gene"] == "PARVB":
            print("\n=== PARVB W37R (AM-missed effector, bonus check) ===")
            print(f"  ESM1b={r.get('esm1b_llr','')} popEVE={r.get('popeve_score','')}"
                  f"  LoGoFunc={r['logofunc_call']} ({r['logofunc_note']})"
                  f"  EVE={r['eve_score']}/{r['eve_class']} ({r['eve_note']})")

    write_readme(rows, gof_dist, logofunc_dist, eve_available, args.eve_merged)


def write_readme(rows, gof_dist, logofunc_dist, eve_available, eve_merged):
    n = len(rows)
    lines = [
        "# coding_hardening_v2.tsv -- coding-arm GoF/LoF + EVE upgrade",
        "",
        "Additive Phase-5 Track-3 output. Extends coding_hardening.tsv "
        f"({n} effectors) with a gain-/loss-of-function direction and an "
        "evolution-only EVE vote. NOMINATION / mechanism-class only; never "
        "summed into the convergence atlas.",
        "",
        "## New columns",
        "- logofunc_neutral / logofunc_gof / logofunc_lof: LoGoFunc "
        "(Stein 2023, Genome Med) precomputed hg38 probabilities "
        "(Zenodo 10126185).",
        "- logofunc_call: argmax(neutral, GOF, LOF) = raw LoGoFunc direction.",
        "- logofunc_note: covered / absent_from_logofunc (position or REF/ALT "
        "not a canonical missense in LoGoFunc).",
        "- gof_lof_call: FINAL direction -- experimentally curated anchor "
        "(curated_literature) if present, else LoGoFunc argmax.",
        "- gof_lof_source: curated_literature / LoGoFunc / none.",
        "- eve_score / eve_class: EVE (Frazer 2021, Nature) evolution-only "
        "severity + Benign/Uncertain/Pathogenic; eve_note = coverage.",
        "",
        "## Held-out validation of the GoF-vs-LoF axis (src/81)",
        "The LoGoFunc GoF-vs-LoF discriminator was benchmarked on LoGoFunc's "
        "OWN official withheld test split (n=2,831; 152 GoF vs 1,340 LoF), "
        "scored by the shipped 27-model ensemble -- NOT on our 42 effectors: "
        "GoF-vs-LoF auROC = **0.905**, 95% CI [0.879, 0.930] (2,000x "
        "bootstrap), permutation-null 0.501, p=5e-4. Lower CI 0.879 >> the "
        "0.70 pre-registered bar -> PASS. CAVEAT: this is VARIANT-level "
        "held-out, NOT gene-family-disjoint (released feature files carry no "
        "gene ID; gene-level features can inflate), so 0.905 is an UPPER BOUND "
        "on strict family-held-out performance; 3-class GoF precision is only "
        "0.49 (GoF rare + over-called among Neutrals). Full numbers: "
        "SF/gof_lof_benchmark.json.",
        "INTERPRETATION RULE: a POSITIVE LoGoFunc GoF/LoF assignment "
        "(logofunc_call in {GOF,LOF}) is a held-out-supported CALL; a Neutral "
        "call on an established mode-of-action gene should DEFER to the "
        "curated_literature flag (see PNPLA3/TM6SF2/APOE below).",
        "",
        "## Why LoGoFunc (not PreMode)",
        "PreMode (Shen 2025) is mode-of-action-native but precomputed for only "
        "10 genes (PNPLA3 not among them); per-gene inference needs AF2 "
        "structure + MSA + GPU. LoGoFunc is genome-wide/turnkey and gives a "
        "GoF/LoF call for every canonical missense -> chosen for coverage.",
        "",
        "## Raw LoGoFunc argmax distribution (all effectors)",
    ]
    for k, v in sorted(logofunc_dist.items(), key=lambda x: -x[1]):
        lines.append(f"- {k}: {v}")
    lines += ["", "## FINAL gof_lof_call distribution (curated-preferred)"]
    for k, v in sorted(gof_dist.items(), key=lambda x: -x[1]):
        lines.append(f"- {k}: {v}")
    lines += [
        "",
        "## PNPLA3 I148M -- the neomorph test",
    ]
    for r in rows:
        if r["gene"] == "PNPLA3" and r["protein_variant"] == "I148M":
            lines.append(
                f"LoGoFunc call = **{r['logofunc_call']}** "
                f"(neutral {r['logofunc_neutral']}, GOF {r['logofunc_gof']}, "
                f"LOF {r['logofunc_lof']}). Final gof_lof_call = "
                f"{r['gof_lof_call']} ({r['gof_lof_source']}). "
                "HONEST: LoGoFunc does NOT recover I148M as GoF (argmax "
                "Neutral; among pathogenic classes it leans LoF). The "
                "curated-literature anchor (retained-on-lipid-droplet "
                "neomorph; BasuRay 2019) remains the correct GoF flag -- "
                "even a GoF/LoF-specialized genome-wide model misses this "
                "neomorph, which is exactly why the curated anchor is kept."
            )
    if not eve_available:
        lines += [
            "",
            "## EVE status",
            f"Merged EVE VCF not present at {eve_merged} at build time "
            "(evemodel.org backend returned 502/504). eve_* columns are blank; "
            "re-run this script (LoGoFunc columns are stable) once EVE lands.",
        ]
    with open(OUT_README, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[out] {OUT_README}")


if __name__ == "__main__":
    main()
