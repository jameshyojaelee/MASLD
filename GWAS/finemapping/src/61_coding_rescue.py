#!/usr/bin/env python3
"""
61_coding_rescue.py  --  Stage-1 "coding rescue" seqfunc annotation.

Annotates the CODING + SPLICE credible-set variants with *precomputed*
deleteriousness so that coding effectors that are invisible to the eQTL/COLOC
layer (e.g. PNPLA3 I148M, TM6SF2 E167K, HSD17B13 splice) are still scored.
This is the non-circular patch that closes the eQTL-blind coding gap: the
scores come from external protein/splice models, NOT from the GWAS/eQTL data
that already defined the credible sets.

Two evidence channels:
  * AlphaMissense (Cheng et al. 2023, Science) -- am_pathogenicity + am_class
    for MISSENSE substitutions, tabix-queried from the precomputed hg38 table.
  * SpliceAI (Jaganathan et al. 2019, Cell) -- delta scores DS_AG/AL/DG/DL
    (acceptor/donor gain/loss) + spliceai_max_ds, computed on the splice /
    splice-region variants with the SpliceAI deep model.

Pipeline (all additive; no existing script touched):
  1. Load credible_set_variant_consequences.csv, keep class in
     {coding_protein_altering, splice_region, coding_synonymous}  (126 vars).
  2. LiftOver hg19 -> hg38 with the repo-standard rtracklayer::liftOver
     (same chain + template as src/55c/58b). If Deliverable-1's substrate
     (results/seqfunc/variant_substrate_hg38.tsv) exists it is used instead.
  3. AlphaMissense join: pysam tabix query per hg38 position, allele-normalise
     to the genomic REF (handles the 5 orientation-flipped duplicate entries),
     assign the missense score of the physical substitution.
  4. SpliceAI: build a de-duplicated hg38 VCF, run the SpliceAI CLI from the
     throwaway env, parse DS_AG/AL/DG/DL and the max.
  5. Emit results/seqfunc/coding_rescue.tsv.

RUN (rnaseq env supplies python+pysam+tabix AND R+rtracklayer):
  micromamba run -n rnaseq python GWAS/finemapping/src/61_coding_rescue.py

Env knobs (all optional):
  MASLD_PROJECT_ROOT   project root override
  SPLICEAI_ENV         path to a python env exposing the `spliceai` CLI
                       (default: <root>/.envs/spliceai_throwaway)
  SKIP_SPLICEAI=1      skip the SpliceAI channel (AlphaMissense only)
"""
import os
import sys
import csv
import shutil
import subprocess
import tempfile

import pysam

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
CS_CSV = os.path.join(ROOT, "GWAS/finemapping/results/credible_set_variant_consequences.csv")
CHAIN = os.path.join(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain")
AM_TSV = os.path.join(ROOT, "data/external/alphamissense/AlphaMissense_hg38.tsv.gz")
GENOME_FA = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

SEQFUNC = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
SUBSTRATE = os.path.join(SEQFUNC, "variant_substrate_hg38.tsv")  # Deliverable 1 (optional)
HG38_TSV = os.path.join(SEQFUNC, "coding_splice_hg38.tsv")       # our own liftOver product
OUT_TSV = os.path.join(SEQFUNC, "coding_rescue.tsv")

SPLICEAI_ENV = os.environ.get("SPLICEAI_ENV", os.path.join(ROOT, ".envs/spliceai_throwaway"))
CODING_SPLICE_CLASSES = {"coding_protein_altering", "splice_region", "coding_synonymous"}

# Resolve Rscript: PATH first, then the rnaseq env (supplies rtracklayer).
RSCRIPT = (shutil.which("Rscript")
           or "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript")
# Set USE_SUBSTRATE=1 to consume Deliverable-1's variant_substrate_hg38.tsv
# instead of our own (validated, decoupled) rtracklayer liftOver.
USE_SUBSTRATE = os.environ.get("USE_SUBSTRATE") == "1"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def log(msg):
    print(f"[61_coding_rescue] {msg}", flush=True)


# ----------------------------------------------------------------------------
# 1. Load coding + splice credible-set variants
# ----------------------------------------------------------------------------
def load_coding_splice():
    rows = []
    with open(CS_CSV) as f:
        for r in csv.DictReader(f):
            if r["class"] in CODING_SPLICE_CLASSES:
                rows.append(
                    {
                        "variant_id": r["variant_id"],          # hg19 chr:pos:ref:alt
                        "chr": r["chr"],
                        "pos": int(r["pos"]),                    # hg19
                        "ref": r["ref"],
                        "alt": r["alt"],
                        "max_pip": r["max_pip"],
                        "gene": r["SYMBOL"],
                        "consequence": r["Consequence"],
                        "var_class": r["class"],
                        "is_indel": r["is_indel"],
                    }
                )
    log(f"loaded {len(rows)} coding+splice credible-set variants")
    return rows


# ----------------------------------------------------------------------------
# 2. LiftOver hg19 -> hg38  (repo-standard rtracklayer, same chain as 55c/58b)
# ----------------------------------------------------------------------------
def liftover_hg38(rows):
    # Optional fast path: reuse Deliverable-1 substrate (schema:
    # variant_id_hg19, chr[=hg38, chr-prefixed], pos_hg38, ...). Off by default
    # so this deliverable stays decoupled from that in-flight artifact.
    if USE_SUBSTRATE and os.path.exists(SUBSTRATE):
        log(f"USE_SUBSTRATE=1 -> reading {SUBSTRATE}")
        sub = {}
        with open(SUBSTRATE) as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for r in rdr:
                key = r.get("variant_id_hg19") or r.get("variant_id") or r.get("variant")
                ch = r.get("chr") or r.get("chr_hg38") or r.get("chrom_hg38")
                po = r.get("pos_hg38") or r.get("hg38_pos")
                if key and ch and po:
                    try:
                        sub[key] = (str(ch).replace("chr", ""), int(po))
                    except ValueError:
                        pass
        hit = sum(1 for r in rows if r["variant_id"] in sub)
        if hit == len(rows):
            for r in rows:
                r["chr_hg38"], r["pos_hg38"] = sub[r["variant_id"]]
            log(f"substrate provided hg38 for all {hit} variants")
            return rows
        log(f"substrate covers only {hit}/{len(rows)} -> falling back to own liftOver")

    # Own liftOver via rtracklayer (rnaseq env provides R + rtracklayer).
    in_tsv = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    for i, r in enumerate(rows):
        in_tsv.write(f"{i}\t{r['chr']}\t{r['pos']}\n")
    in_tsv.close()
    r_script = tempfile.NamedTemporaryFile("w", suffix=".R", delete=False)
    r_script.write(
        f"""
suppressMessages({{library(rtracklayer); library(GenomicRanges)}})
d <- read.table("{in_tsv.name}", sep="\\t", col.names=c("idx","chr","pos"),
                colClasses=c("integer","character","integer"))
chain <- import.chain("{CHAIN}")
gr <- GRanges(paste0("chr", d$chr), IRanges(d$pos, width=1)); gr$idx <- d$idx
lifted <- liftOver(gr, chain)
n <- lengths(lifted)
out <- data.frame(idx=integer(0), chr_hg38=character(0), pos_hg38=integer(0))
for (i in seq_along(lifted)) {{
  if (length(lifted[[i]]) == 1) {{
    out <- rbind(out, data.frame(idx=d$idx[i],
      chr_hg38=sub("^chr","",as.character(seqnames(lifted[[i]]))),
      pos_hg38=start(lifted[[i]])))
  }}
}}
write.table(out, "{HG38_TSV}", sep="\\t", quote=FALSE, row.names=FALSE)
cat("lifted", nrow(out), "of", nrow(d), "\\n")
"""
    )
    r_script.close()
    log(f"running rtracklayer liftOver hg19->hg38 via {RSCRIPT} ...")
    subprocess.run([RSCRIPT, r_script.name], check=True)

    lifted = {}
    with open(HG38_TSV) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for rr in rdr:
            lifted[int(rr["idx"])] = (rr["chr_hg38"], int(rr["pos_hg38"]))
    for i, r in enumerate(rows):
        if i in lifted:
            r["chr_hg38"], r["pos_hg38"] = lifted[i]
        else:
            r["chr_hg38"], r["pos_hg38"] = None, None
    n_lift = sum(1 for r in rows if r["pos_hg38"] is not None)
    log(f"liftOver: {n_lift}/{len(rows)} lifted 1:1")
    os.unlink(in_tsv.name)
    os.unlink(r_script.name)
    return rows


# ----------------------------------------------------------------------------
# 3. AlphaMissense join  (allele-normalised to genomic REF)
# ----------------------------------------------------------------------------
def annotate_alphamissense(rows):
    if not os.path.exists(AM_TSV) or not os.path.exists(AM_TSV + ".tbi"):
        log(f"FATAL: AlphaMissense table or index missing at {AM_TSV}")
        sys.exit(2)
    am = pysam.TabixFile(AM_TSV)
    am_contigs = set(am.contigs)

    for r in rows:
        r.update(dict(am_pathogenicity="", am_class="", am_protein_variant="",
                      am_transcript="", am_note=""))
        if r["pos_hg38"] is None:
            r["am_note"] = "no_hg38"
            continue
        contig = "chr" + str(r["chr_hg38"]) if not str(r["chr_hg38"]).startswith("chr") else str(r["chr_hg38"])
        if contig not in am_contigs:
            r["am_note"] = "contig_absent"
            continue
        pos = r["pos_hg38"]
        recs = []
        for line in am.fetch(contig, pos - 1, pos):
            f = line.split("\t")
            # CHROM POS REF ALT genome uniprot transcript protein_var am_path am_class
            if int(f[1]) != pos:
                continue
            recs.append(f)
        if not recs:
            r["am_note"] = "not_missense_in_AM"   # synonymous / splice / no coverage
            continue
        gref = recs[0][2]  # genomic REF (+ strand), identical across rows at a pos
        a1, a2 = r["ref"], r["alt"]
        # physical alt = allele that is NOT the genomic reference
        if a1 == gref:
            sub_alt, orient = a2, "fwd"
        elif a2 == gref:
            sub_alt, orient = a1, "flip"
        elif COMP.get(a1) == gref:
            sub_alt, orient = COMP.get(a2), "fwd_strandflip"
        elif COMP.get(a2) == gref:
            sub_alt, orient = COMP.get(a1), "flip_strandflip"
        else:
            r["am_note"] = f"allele_mismatch(gref={gref};id={a1}/{a2})"
            continue
        matched = [f for f in recs if f[3] == sub_alt]
        if not matched:
            r["am_note"] = f"alt_not_scored(sub_alt={sub_alt})"
            continue
        # pick most-deleterious transcript row; record how many transcripts
        best = max(matched, key=lambda f: float(f[8]))
        r["am_pathogenicity"] = best[8]
        r["am_class"] = best[9]
        r["am_protein_variant"] = best[7]
        r["am_transcript"] = best[6]
        note = [] if orient == "fwd" else [f"orient={orient}"]
        if len(matched) > 1:
            note.append(f"n_tx={len(matched)}")
        if r["var_class"] != "coding_protein_altering":
            note.append("AM_missense_vs_VEP_" + r["var_class"])
        r["am_note"] = ";".join(note)
    am.close()
    n_scored = sum(1 for r in rows if r["am_pathogenicity"] != "")
    log(f"AlphaMissense: scored {n_scored}/{len(rows)} variants")
    return rows


# ----------------------------------------------------------------------------
# 4. SpliceAI  (build de-duplicated hg38 VCF, run CLI, parse DS)
# ----------------------------------------------------------------------------
def run_spliceai(rows):
    for r in rows:
        r.update(dict(spliceai_max_ds="", spliceai_DS_AG="", spliceai_DS_AL="",
                      spliceai_DS_DG="", spliceai_DS_DL="", spliceai_note=""))
    if os.environ.get("SKIP_SPLICEAI") == "1":
        for r in rows:
            r["spliceai_note"] = "skipped"
        log("SpliceAI: SKIP_SPLICEAI=1 -> deferred")
        return rows, "skipped"

    spliceai_bin = os.path.join(SPLICEAI_ENV, "bin", "spliceai")
    if not os.path.exists(spliceai_bin):
        for r in rows:
            r["spliceai_note"] = "env_unavailable"
        log(f"SpliceAI: CLI not found at {spliceai_bin} -> deferred (partial)")
        return rows, "env_unavailable"

    # Build genomic-REF-normalised, de-duplicated VCF over splice-relevant variants.
    # Run SpliceAI on every coding+splice variant that lifted (cheap); DS~0 for
    # deep-exonic sites, informative near splice junctions.
    fa = pysam.FastaFile(GENOME_FA)
    seen = {}          # (chrom,pos,ref,alt) -> vcf_id
    key_for_var = {}   # variant_id -> (chrom,pos,ref,alt)
    vcf_records = []
    for r in rows:
        if r["pos_hg38"] is None:
            continue
        contig = "chr" + str(r["chr_hg38"]) if not str(r["chr_hg38"]).startswith("chr") else str(r["chr_hg38"])
        pos = r["pos_hg38"]
        try:
            gref = fa.fetch(contig, pos - 1, pos).upper()
        except Exception:
            r["spliceai_note"] = "fasta_fetch_fail"
            continue
        a1, a2 = r["ref"], r["alt"]
        if a1 == gref:
            alt = a2
        elif a2 == gref:
            alt = a1
        elif COMP.get(a1) == gref:
            alt = COMP.get(a2)
        elif COMP.get(a2) == gref:
            alt = COMP.get(a1)
        else:
            r["spliceai_note"] = f"ref_mismatch(gref={gref})"
            continue
        if alt is None or alt == gref:
            r["spliceai_note"] = "no_alt"
            continue
        key = (contig, pos, gref, alt)
        key_for_var[r["variant_id"]] = key
        if key not in seen:
            seen[key] = f"{contig}_{pos}_{gref}_{alt}"
            vcf_records.append(key)
    fa.close()

    if not vcf_records:
        log("SpliceAI: no records to score")
        return rows, "no_records"

    vcf_in = os.path.join(SEQFUNC, "spliceai_input.vcf")
    vcf_out = os.path.join(SEQFUNC, "spliceai_output.vcf")
    contigs_sorted = sorted({k[0] for k in vcf_records},
                            key=lambda c: (len(c), c))
    with open(vcf_in, "w") as vf:
        vf.write("##fileformat=VCFv4.2\n")
        for c in contigs_sorted:
            vf.write(f"##contig=<ID={c}>\n")
        vf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (contig, pos, gref, alt) in sorted(vcf_records, key=lambda k: (len(k[0]), k[0], k[1])):
            vf.write(f"{contig}\t{pos}\t{contig}_{pos}_{gref}_{alt}\t{gref}\t{alt}\t.\t.\t.\n")
    log(f"SpliceAI: wrote {len(vcf_records)} records to {vcf_in}")

    cmd = [spliceai_bin, "-I", vcf_in, "-O", vcf_out,
           "-R", GENOME_FA, "-A", "grch38"]
    log("SpliceAI: running " + " ".join(cmd))
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", TF_CPP_MIN_LOG_LEVEL="3")
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.exists(vcf_out):
        log(f"SpliceAI FAILED rc={proc.returncode}; stderr tail:\n{proc.stderr[-1500:]}")
        for r in rows:
            if not r["spliceai_note"]:
                r["spliceai_note"] = "spliceai_run_failed"
        return rows, "run_failed"

    # Parse SpliceAI INFO: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    ds_by_key = {}
    with open(vcf_out) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            info = f[7]
            if "SpliceAI=" not in info:
                continue
            best = None
            for field in info.split(";"):
                if not field.startswith("SpliceAI="):
                    continue
                for ann in field[len("SpliceAI="):].split(","):
                    p = ann.split("|")
                    if len(p) < 6:
                        continue
                    try:
                        ds = [float(p[2]), float(p[3]), float(p[4]), float(p[5])]
                    except ValueError:
                        continue
                    if best is None or max(ds) > max(best[1]):
                        best = (p[1], ds)  # (symbol, [AG,AL,DG,DL])
            if best is not None:
                ds_by_key[(chrom, pos, ref, alt)] = best[1]
    for r in rows:
        key = key_for_var.get(r["variant_id"])
        if key and key in ds_by_key:
            ag, al, dg, dl = ds_by_key[key]
            r["spliceai_DS_AG"] = f"{ag:.4f}"
            r["spliceai_DS_AL"] = f"{al:.4f}"
            r["spliceai_DS_DG"] = f"{dg:.4f}"
            r["spliceai_DS_DL"] = f"{dl:.4f}"
            r["spliceai_max_ds"] = f"{max(ag, al, dg, dl):.4f}"
            if not r["spliceai_note"]:
                r["spliceai_note"] = "scored"
        elif key and not r["spliceai_note"]:
            r["spliceai_note"] = "not_in_annotation"
    n_scored = sum(1 for r in rows if r["spliceai_max_ds"] != "")
    log(f"SpliceAI: scored {n_scored}/{len(rows)} variants")
    return rows, "ok"


# ----------------------------------------------------------------------------
# 5. Write output
# ----------------------------------------------------------------------------
def write_output(rows):
    cols = [
        "variant", "variant_hg38", "gene", "var_class", "consequence", "max_pip",
        "am_pathogenicity", "am_class", "am_protein_variant", "am_transcript", "am_note",
        "spliceai_max_ds", "spliceai_DS_AG", "spliceai_DS_AL",
        "spliceai_DS_DG", "spliceai_DS_DL", "spliceai_note",
    ]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            v38 = (f"{r['chr_hg38']}:{r['pos_hg38']}:{r['ref']}:{r['alt']}"
                   if r["pos_hg38"] is not None else "")
            w.writerow([
                r["variant_id"], v38, r["gene"], r["var_class"], r["consequence"], r["max_pip"],
                r["am_pathogenicity"], r["am_class"], r["am_protein_variant"],
                r["am_transcript"], r["am_note"],
                r["spliceai_max_ds"], r["spliceai_DS_AG"], r["spliceai_DS_AL"],
                r["spliceai_DS_DG"], r["spliceai_DS_DL"], r["spliceai_note"],
            ])
    log(f"wrote {OUT_TSV} ({len(rows)} rows)")


def main():
    os.makedirs(SEQFUNC, exist_ok=True)
    rows = load_coding_splice()
    rows = liftover_hg38(rows)
    rows = annotate_alphamissense(rows)
    rows, sp_status = run_spliceai(rows)
    write_output(rows)

    # ---- verification summary --------------------------------------------
    pa = [r for r in rows if r["var_class"] == "coding_protein_altering"]
    pa_genes = sorted({r["gene"] for r in pa})
    pa_genes_scored = sorted({r["gene"] for r in pa if r["am_pathogenicity"] != ""})
    log("=" * 60)
    log(f"coding_protein_altering genes: {len(pa_genes)}; "
        f"with >=1 AM-scored variant: {len(pa_genes_scored)}")
    missing = set(pa_genes) - set(pa_genes_scored)
    if missing:
        log(f"protein-altering genes with NO AM score: {sorted(missing)}")
    pnpla3 = [r for r in rows if r["gene"] == "PNPLA3"
              and r["var_class"] == "coding_protein_altering"
              and r["pos_hg38"] == 43928847]
    if pnpla3:
        p = pnpla3[0]
        log(f"PNPLA3 I148M (rs738409, chr22:43928847 C>G): "
            f"am_pathogenicity={p['am_pathogenicity']} am_class={p['am_class']} "
            f"protein={p['am_protein_variant']} tx={p['am_transcript']}")
    n_am = sum(1 for r in rows if r["am_pathogenicity"] != "")
    n_sp = sum(1 for r in rows if r["spliceai_max_ds"] != "")
    log(f"TOTAL: {len(rows)} variants | AM scored {n_am} | SpliceAI scored {n_sp} "
        f"| SpliceAI status: {sp_status}")


if __name__ == "__main__":
    main()
