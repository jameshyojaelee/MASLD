#!/usr/bin/env python3
"""
66_coding_hardening.py  --  Axis-1 CODING-ARM HARDENING for the seqfunc layer.

Closes the PARVB gap that AlphaMissense (AM) missed and puts all 34 coding
effectors on ONE cross-gene-calibrated severity scale by adding two protein
language-model channels that AM does NOT provide:

  * ESM1b  (Brandes et al. 2023, Nat Genet) -- masked-marginal LLR of the
    substitution, log P(mut|context) - log P(wt|context). Same 650M model /
    same scale for every gene, so ESM LLR is comparable across genes in a way
    AM's per-gene-uncalibrated pathogenicity is not.
  * popEVE (Orenbuch/Marks 2025, medRxiv 2023.11.27.23299062) -- a
    PROTEOME-WIDE, human-population-calibrated severity that places every
    missense on ONE cross-gene spectrum. This is the key upgrade over AM:
    AM pathogenicity is NOT cross-gene calibrated, popEVE is.

Design (STRICTLY ADDITIVE; no existing 46d/78/27a/60-65 script is touched):
  mode=prep     Build the effector variant table from coding_rescue.tsv
                (34 coding_protein_altering missense genes) + PARVB, resolve a
                protein sequence + WT/pos/MUT residue for each (gencode v49
                translations, with per-gene canonical fallback for version-
                mismatched AM transcripts).  ->  coding_hardening_prep.tsv
  mode=esm      Run ESM1b masked-marginal LLR (GPU, fair-esm, offline weights).
                ->  coding_hardening_esm.tsv
  mode=popeve   Coordinate lookup of the popEVE GRCh38 VCF (tabix).
                ->  coding_hardening_popeve.tsv
  mode=assemble Join AM(from coding_rescue) + ESM + popEVE + curated GoF/LoF,
                rank cross-gene by popEVE, flag AM-benign/ESM-or-popEVE-damaging
                disagreements.  ->  coding_hardening.tsv  (final deliverable)
  mode=all      prep -> popeve -> assemble  (ESM read from esm.tsv if present).

MERE NOMINATION / MECHANISM-CLASS ONLY. These columns are protein-model
severity for the eQTL-blind coding effectors; they are NEVER wired into the
convergence atlas or a scored channel (would be circular with COLOC / epi).

RUN
  # 1. prep + popeve (CPU, fast) in rnaseq env (pysam+tabix+Bio)
  micromamba run -n rnaseq python GWAS/finemapping/src/66_coding_hardening.py --mode prep
  micromamba run -n rnaseq python GWAS/finemapping/src/66_coding_hardening.py --mode popeve
  # 2. ESM1b on GPU
  sbatch GWAS/finemapping/src/66_esm.sbatch
  # 3. assemble
  micromamba run -n rnaseq python GWAS/finemapping/src/66_coding_hardening.py --mode assemble
"""
import os
import re
import csv
import gzip
import glob
import argparse
import subprocess

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SEQFUNC = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
CODING_RESCUE = os.path.join(SEQFUNC, "coding_rescue.tsv")
SUBSTRATE = os.path.join(SEQFUNC, "variant_substrate_hg38.tsv")
GENCODE_PEP = os.path.join(ROOT, "data/external/gencode_v49_pep/gencode.v49.pc_translations.fa.gz")
GENCODE_GTF = "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
GENOME_FA = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
POPEVE_VCF = os.path.join(ROOT, "data/external/popeve/grch38_popEVE_ukbb_20250715.vcf.gz")

PREP_TSV = os.path.join(SEQFUNC, "coding_hardening_prep.tsv")
ESM_TSV = os.path.join(SEQFUNC, "coding_hardening_esm.tsv")
POPEVE_TSV = os.path.join(SEQFUNC, "coding_hardening_popeve.tsv")
OUT_TSV = os.path.join(SEQFUNC, "coding_hardening.tsv")

ESM_TORCH_HOME = os.path.join(ROOT, "data/external/esm_weights/torch_home")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ---------------------------------------------------------------------------
# PARVB special case: AM is MANE-restricted (ENST00000338758 / NP_...), and the
# credible-set missense sits UPSTREAM of the MANE TSS -> non-coding in MANE, so
# AM returns nothing. It IS a genuine missense on the alternative-N-terminal
# isoform ENST00000406477.7 (RefSeq NP_001003828.1): chr22:43999571 T>C = W37R
# (ref codon TGG -> alt codon CGG). We route it through ESM/popEVE on that
# isoform. (Diagnosis verified: gencode ENSP residue 37 = W; popEVE VCF row at
# 22:43999571 T>C reports gene=PARVB mutant=W37R.)
PARVB_ROW = dict(
    variant="22:44395451:T:C",
    variant_hg38="22:43999571:T:C",
    chr="22", pos_hg38="43999571", ref="T", alt="C",
    gene="PARVB", transcript="ENST00000406477.7", protein_id="NP_001003828.1",
    protein_variant="W37R", wt="W", aa_pos="37", mut="R",
    seq_source="gencode_v49_alt_isoform",
    am_pathogenicity="", am_class="",
    prep_note="AM_missed=MANE_upstream;missense_on_alt_isoform_ENST00000406477",
)

# ---------------------------------------------------------------------------
# Curated GoF / LoF flags (LoGoFunc bulk download unavailable: repo 404,
# Zenodo record 410, only an interactive Shiny app). Literature-curated for the
# established effectors; clearly tagged source=curated_literature, NOT computed.
CURATED_GOF_LOF = {
    # gene -> (protein_variant, flag, note)
    "PNPLA3": ("I148M", "GoF_neomorph",
               "rs738409 I148M: retained-on-lipid-droplet neomorph / loss of "
               "hydrolase + gain of aberrant LD sequestration (Smagris 2015; "
               "BasuRay 2019). AM=0.32 benign-ish underit; functional GoF."),
    "TM6SF2": ("E167K", "LoF",
               "rs58542926 E167K: destabilising LoF, reduced VLDL secretion "
               "(Kozlitina 2014; Smagris 2016)."),
    "HSD17B13": ("", "LoF_protective",
                 "rs72613567 splice: LoF, protective (Abul-Husn 2018) -- "
                 "splice not missense; not in this missense arm."),
    "GCKR": ("P446L", "functional_common",
             "rs1260326 P446L: common functional variant, altered GCKR-F6P "
             "regulation (Beer 2009)."),
    "MARC1": ("A165T", "LoF_protective",
              "rs2642438 A165T: protective missense (Emdin 2020)."),
    "APOE": ("", "isoform_functional",
             "APOE isoform-defining residue (C130R / R176C axis)."),
}


def log(m):
    print(f"[66_coding_hardening] {m}", flush=True)


# ===========================================================================
# gencode protein index
# ===========================================================================
def load_gencode_pep():
    """Return (by_enst_versioned, by_enst_base, by_gene[gene]->[(enst,seq)])."""
    pepv, pepb, byg = {}, {}, {}

    def flush(hdr, seq):
        if not hdr:
            return
        parts = hdr[1:].split("|")
        ensts = [p for p in parts if p.startswith("ENST")]
        gname = parts[6] if len(parts) > 6 else parts[-1]
        for e in ensts:
            pepv[e] = seq
            pepb.setdefault(e.split(".")[0], seq)
        if ensts:
            byg.setdefault(gname, []).append((ensts[0], seq))

    with gzip.open(GENCODE_PEP, "rt") as f:
        hdr, seq = None, []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                flush(hdr, "".join(seq))
                hdr, seq = line, []
            else:
                seq.append(line)
        flush(hdr, "".join(seq))
    return pepv, pepb, byg


# ===========================================================================
# mode=prep
# ===========================================================================
def mode_prep():
    pepv, pepb, byg = load_gencode_pep()
    log(f"gencode v49: {len(pepb)} base-ENST protein sequences")

    rows = []
    with open(CODING_RESCUE) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["var_class"] == "coding_protein_altering" and "missense" in r["consequence"]:
                # PARVB is handled explicitly below (its coding_rescue row is the
                # blank AM-missed entry); skip it here to avoid a duplicate.
                if r["gene"] == "PARVB":
                    continue
                rows.append(r)
    log(f"coding_protein_altering missense rows in coding_rescue (excl PARVB): {len(rows)} "
        f"({len({r['gene'] for r in rows})} genes)")

    out = []
    for r in rows:
        gene = r["gene"]
        pv = r["am_protein_variant"]
        tx = r["am_transcript"]
        rec = dict(
            variant=r["variant"], variant_hg38=r["variant_hg38"],
            gene=gene, transcript=tx, protein_id="", protein_variant=pv,
            wt="", aa_pos="", mut="", seq_source="", prep_note="",
            am_pathogenicity=r["am_pathogenicity"], am_class=r["am_class"],
        )
        # parse variant_hg38 "chr:pos:ref:alt"
        p = r["variant_hg38"].split(":")
        rec["chr"], rec["pos_hg38"], rec["ref"], rec["alt"] = p[0], p[1], p[2], p[3]

        m = re.match(r"^([A-Z])(\d+)([A-Z])$", pv or "")
        if not m:
            rec["prep_note"] = "unparseable_protein_variant"
            out.append(rec)
            continue
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        rec["wt"], rec["aa_pos"], rec["mut"] = wt, str(pos), mut

        # resolve a protein sequence whose residue `pos` == wt
        seq, src = None, ""
        s = pepv.get(tx) or pepb.get(tx.split(".")[0])
        if s and pos - 1 < len(s) and s[pos - 1] == wt:
            seq, src = s, "gencode_v49_am_transcript"
        else:
            # per-gene fallback: any gencode transcript of this gene with wt@pos
            for enst, cand in byg.get(gene, []):
                if pos - 1 < len(cand) and cand[pos - 1] == wt:
                    seq, src, tx = cand, "gencode_v49_gene_fallback", enst
                    rec["transcript"] = enst
                    break
        if seq is None:
            rec["seq_source"] = "UNRESOLVED"
            rec["prep_note"] = ("am_transcript_retired_in_gencode_v49;"
                                "esm_not_computable(popEVE_covers_by_coord)")
        else:
            rec["seq_source"] = src
        rec["_seq"] = seq or ""
        out.append(rec)

    # add PARVB (resolve its isoform sequence for ESM)
    parvb = dict(PARVB_ROW)
    s = pepv.get(parvb["transcript"])
    if s and int(parvb["aa_pos"]) - 1 < len(s) and s[int(parvb["aa_pos"]) - 1] == parvb["wt"]:
        parvb["_seq"] = s
    else:
        parvb["_seq"] = ""
        parvb["seq_source"] = "UNRESOLVED"
    out.append(parvb)

    cols = ["gene", "variant", "variant_hg38", "chr", "pos_hg38", "ref", "alt",
            "transcript", "protein_id", "protein_variant", "wt", "aa_pos", "mut",
            "seq_source", "am_pathogenicity", "am_class", "prep_note", "seq"]
    with open(PREP_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in out:
            w.writerow([r.get(c if c != "seq" else "_seq", "") for c in cols])
    n_seq = sum(1 for r in out if r.get("_seq"))
    log(f"wrote {PREP_TSV}: {len(out)} variants, {n_seq} with a protein sequence "
        f"(ESM-computable), {len(out)-n_seq} unresolved")
    return out


# ===========================================================================
# mode=esm   (GPU, fair-esm masked-marginal LLR)
# ===========================================================================
def mode_esm():
    import torch
    os.environ.setdefault("TORCH_HOME", ESM_TORCH_HOME)
    import esm

    rows = []
    with open(PREP_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    todo = [r for r in rows if r["seq"] and r["wt"] and r["aa_pos"] and r["mut"]]
    log(f"ESM1b: {len(todo)}/{len(rows)} variants have a sequence")

    log("loading esm1b_t33_650M_UR50S (offline from TORCH_HOME) ...")
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    bc = alphabet.get_batch_converter()
    model.eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev)
    log(f"device={dev}")
    MAXRES = 1022  # ESM1b context = 1024 tokens incl BOS/EOS

    results = {}
    for r in todo:
        seq = r["seq"]
        pos = int(r["aa_pos"])          # 1-based in full protein
        wt, mut = r["wt"], r["mut"]
        # sliding window centred on the variant for long proteins
        if len(seq) <= MAXRES:
            wseq, lpos = seq, pos
        else:
            start = max(0, (pos - 1) - MAXRES // 2)
            start = min(start, len(seq) - MAXRES)
            wseq = seq[start:start + MAXRES]
            lpos = (pos - 1) - start + 1
        assert wseq[lpos - 1] == wt, f"{r['gene']} window WT mismatch"

        # mask the variant position -> masked-marginal
        masked = list(wseq)
        _, _, toks = bc([("v", "".join(masked))])
        toks = toks.to(dev)
        toks[0, lpos] = alphabet.mask_idx  # +1 for BOS token
        with torch.no_grad():
            logits = model(toks)["logits"]
        lp = torch.log_softmax(logits[0, lpos], dim=-1)
        llr = (lp[alphabet.get_idx(mut)] - lp[alphabet.get_idx(wt)]).item()
        results[r["variant_hg38"] + "|" + r["gene"]] = round(llr, 4)
        log(f"  {r['gene']:9s} {wt}{pos}{mut}  ESM1b_LLR={llr:+.3f}"
            f"{'  [windowed]' if len(seq) > MAXRES else ''}")

    with open(ESM_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["variant_hg38", "gene", "esm1b_llr"])
        for r in rows:
            k = r["variant_hg38"] + "|" + r["gene"]
            if k in results:
                w.writerow([r["variant_hg38"], r["gene"], results[k]])
    log(f"wrote {ESM_TSV}: {len(results)} ESM1b LLR scores")


# ===========================================================================
# mode=popeve   (tabix coordinate lookup of the popEVE GRCh38 VCF)
# ===========================================================================
def _parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def mode_popeve():
    import pysam
    if not os.path.exists(POPEVE_VCF + ".tbi"):
        pysam.tabix_index(POPEVE_VCF, preset="vcf", force=True)
    tb = pysam.TabixFile(POPEVE_VCF)

    rows = []
    with open(PREP_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)

    out = []
    for r in rows:
        chrom, pos = r["chr"].replace("chr", ""), int(r["pos_hg38"])
        ref, alt = r["ref"], r["alt"]
        rec = dict(variant_hg38=r["variant_hg38"], gene=r["gene"],
                   popeve_score="", popeve_mutant="", popeve_protein="",
                   popeve_eve="", popeve_esm1v="", popeve_note="")
        hits = []
        try:
            for line in tb.fetch(chrom, pos - 1, pos):
                f7 = line.split("\t")
                if int(f7[1]) != pos:
                    continue
                hits.append(f7)
        except (ValueError, IndexError):
            rec["popeve_note"] = "contig_absent"
        # Resolve the physical substitution against the popEVE genomic REF, then
        # match the enumerated row with that exact ALT (mirrors script 61's
        # allele-normalisation). Credible-set ref/alt may be swapped and/or
        # strand-flipped relative to the + strand genome; anchoring on the VCF's
        # own REF avoids grabbing a reverse-complement row (e.g. minus-strand
        # TMC4: alleles C/T at a genome-T site -> physical T>C = E17G, not T>A).
        matched = None
        if hits:
            vref = hits[0][3]  # genomic REF, identical across rows at this pos
            if ref == vref:
                phys_alt = alt
            elif alt == vref:
                phys_alt = ref
            elif COMP.get(ref) == vref:
                phys_alt = COMP.get(alt)
            elif COMP.get(alt) == vref:
                phys_alt = COMP.get(ref)
            else:
                phys_alt = None
                rec["popeve_note"] = f"allele_mismatch(vref={vref};id={ref}/{alt})"
            if phys_alt is not None:
                for f7 in hits:
                    if f7[3] == vref and f7[4] == phys_alt:
                        matched = f7; break
        if matched is not None:
            info = _parse_info(matched[7])
            rec["popeve_score"] = info.get("popEVE", "")
            rec["popeve_mutant"] = info.get("mutant", "")
            rec["popeve_protein"] = info.get("protein", "")
            rec["popeve_eve"] = info.get("EVE", "")
            rec["popeve_esm1v"] = info.get("ESM1v", "")
            rec["popeve_note"] = "matched"
            if info.get("gene") and info["gene"] != r["gene"]:
                rec["popeve_note"] = f"matched;gene_annot={info['gene']}"
        elif not rec["popeve_note"]:
            rec["popeve_note"] = "no_popEVE_at_coord" if not hits else "alt_not_in_popEVE"
        out.append(rec)
    tb.close()

    cols = ["variant_hg38", "gene", "popeve_score", "popeve_mutant",
            "popeve_protein", "popeve_eve", "popeve_esm1v", "popeve_note"]
    with open(POPEVE_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in out:
            w.writerow([r[c] for c in cols])
    n = sum(1 for r in out if r["popeve_score"])
    log(f"wrote {POPEVE_TSV}: {n}/{len(out)} variants matched in popEVE VCF")
    return out


# ===========================================================================
# mode=assemble
# ===========================================================================
def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mode_assemble():
    prep = {}
    order = []
    with open(PREP_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = r["variant_hg38"] + "|" + r["gene"]
            prep[k] = r
            order.append(k)

    esm = {}
    if os.path.exists(ESM_TSV):
        with open(ESM_TSV) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                esm[r["variant_hg38"] + "|" + r["gene"]] = r["esm1b_llr"]
    else:
        log(f"WARNING: {ESM_TSV} absent -> esm1b_llr blank (run --mode esm first)")

    pop = {}
    if os.path.exists(POPEVE_TSV):
        with open(POPEVE_TSV) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                pop[r["variant_hg38"] + "|" + r["gene"]] = r

    recs = []
    for k in order:
        p = prep[k]
        gene = p["gene"]
        am = _f(p["am_pathogenicity"])
        esm_llr = _f(esm.get(k, ""))
        popr = pop.get(k, {})
        popeve = _f(popr.get("popeve_score", ""))

        # QC: AA identity of the substitution should agree between the AM-derived
        # residue and the popEVE RefSeq annotation (position numbering may differ
        # by transcript, but WT/MUT amino acids must match; a mismatch would mean
        # a mis-mapped residue in a gene-fallback transcript).
        aa_qc = ""
        pm = re.match(r"^([A-Z])(\d+)([A-Z])$", popr.get("popeve_mutant", "") or "")
        if pm and p["wt"] and p["mut"]:
            if (pm.group(1), pm.group(3)) != (p["wt"], p["mut"]):
                aa_qc = (f"AA_mismatch_vs_popEVE(am={p['wt']}>{p['mut']},"
                         f"popEVE={pm.group(1)}>{pm.group(3)})")

        # curated GoF/LoF (literature; LoGoFunc bulk unavailable)
        gof = ""
        cur = CURATED_GOF_LOF.get(gene)
        if cur and (cur[0] == "" or cur[0] == p["protein_variant"]):
            gof = f"{cur[1]} [curated_literature]"

        # AM class flag
        am_call = ""
        if am is not None:
            am_call = "AM_lpath" if am >= 0.564 else ("AM_lbenign" if am <= 0.34 else "AM_ambiguous")

        # ESM / popEVE damaging calls (directional: more negative = more damaging)
        # ESM1b clinical operating point ~ LLR < -7.5 (Brandes 2023); use tiers.
        esm_call = ""
        if esm_llr is not None:
            esm_call = ("ESM_damaging" if esm_llr <= -7.5 else
                        "ESM_moderate" if esm_llr <= -4.0 else "ESM_tolerant")
        # popEVE proteome-wide: more negative = more pathogenic (ranked below).
        pop_call = ""
        if popeve is not None:
            pop_call = ("popEVE_damaging" if popeve <= -5.0 else
                        "popEVE_moderate" if popeve <= -3.0 else "popEVE_tolerant")

        # disagreement: AM-benign/ambiguous/absent but a protein-LM flags damage.
        # These credible-set effectors are COMMON variants, so popEVE (proteome-
        # wide) is compressed into a moderate band and rarely hits its damaging
        # tier; the gene-comparable ESM1b LLR carries the actionable signal, so
        # ESM_moderate (LLR<=-4) counts as an AM-undercall flag.
        disagree = ""
        if am_call in ("AM_lbenign", "AM_ambiguous") or am is None:
            if esm_call in ("ESM_damaging", "ESM_moderate") or pop_call == "popEVE_damaging":
                disagree = "AM_undercall_vs_LM"

        notes = []
        if p["prep_note"]:
            notes.append(p["prep_note"])
        if p["gene"] == "PARVB":
            notes.append("AM_gap_closed_by_ESM+popEVE")
        if disagree:
            notes.append(disagree)
        if gof:
            notes.append("GoF/LoF=" + gof)
        pn = popr.get("popeve_note", "")
        if pn and pn != "matched":
            notes.append("popEVE:" + pn)
        if aa_qc:
            notes.append(aa_qc)

        recs.append(dict(
            gene=gene, variant=p["variant"], variant_hg38=p["variant_hg38"],
            protein_variant=p["protein_variant"], transcript=p["transcript"],
            esm1b_llr=("" if esm_llr is None else f"{esm_llr:.4f}"),
            popeve_score=("" if popeve is None else f"{popeve:.4f}"),
            am_pathogenicity=p["am_pathogenicity"], am_class=p["am_class"],
            gof_lof_flag=gof,
            esm_call=esm_call, popeve_call=pop_call, am_call=am_call,
            consensus_note=";".join(notes),
            seq_source=p["seq_source"],
            _popeve=popeve, _esm=esm_llr, _am=am,
        ))

    # cross-gene severity ranking by popEVE (the cross-gene-calibrated axis);
    # ascending = most damaging first. Variants without popEVE ranked last.
    ranked = sorted(recs, key=lambda r: (r["_popeve"] is None, r["_popeve"] if r["_popeve"] is not None else 0))
    for i, r in enumerate(ranked, 1):
        r["popeve_severity_rank"] = i if r["_popeve"] is not None else ""

    cols = ["gene", "variant", "variant_hg38", "protein_variant", "transcript",
            "esm1b_llr", "popeve_score", "popeve_severity_rank",
            "am_pathogenicity", "am_class", "gof_lof_flag",
            "esm_call", "popeve_call", "am_call", "seq_source", "consensus_note"]
    # write in original (gene) order for the file; rank column carries the ranking
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in recs:
            w.writerow([r[c] for c in cols])
    log(f"wrote {OUT_TSV}: {len(recs)} coding effectors")

    # ---- report ----------------------------------------------------------
    log("=" * 70)
    n_am = sum(1 for r in recs if r["_am"] is not None)
    n_esm = sum(1 for r in recs if r["_esm"] is not None)
    n_pop = sum(1 for r in recs if r["_popeve"] is not None)
    log(f"coverage: AM {n_am}/{len(recs)} | ESM1b {n_esm}/{len(recs)} | popEVE {n_pop}/{len(recs)}")
    log("--- cross-gene severity ranking (popEVE ascending = most damaging) ---")
    for r in ranked:
        log(f"  #{str(r['popeve_severity_rank']):>3} {r['gene']:9s} "
            f"{r['protein_variant']:8s} popEVE={str(r['popeve_score']):>8} "
            f"ESM1b={str(r['esm1b_llr']):>9} AM={str(r['am_pathogenicity']):>7}")
    dis = [r for r in recs if "AM_undercall_vs_LM" in r["consensus_note"]]
    log(f"--- AM-benign/ambiguous but LM-damaging (n={len(dis)}) ---")
    for r in dis:
        log(f"  {r['gene']:9s} {r['protein_variant']:8s} "
            f"AM={r['am_pathogenicity']}({r['am_class']}) "
            f"ESM1b={r['esm1b_llr']} popEVE={r['popeve_score']}")
    parvb = [r for r in recs if r["gene"] == "PARVB"]
    if parvb:
        r = parvb[0]
        log("--- PARVB (AM gap) ---")
        log(f"  {r['protein_variant']} on {r['transcript']}: "
            f"ESM1b={r['esm1b_llr']} popEVE={r['popeve_score']} "
            f"rank={r['popeve_severity_rank']}/{n_pop}  note={r['consensus_note']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="all",
                    choices=["prep", "esm", "popeve", "assemble", "all"])
    a = ap.parse_args()
    os.makedirs(SEQFUNC, exist_ok=True)
    if a.mode == "prep":
        mode_prep()
    elif a.mode == "esm":
        mode_esm()
    elif a.mode == "popeve":
        mode_popeve()
    elif a.mode == "assemble":
        mode_assemble()
    elif a.mode == "all":
        mode_prep()
        mode_popeve()
        mode_assemble()


if __name__ == "__main__":
    main()
