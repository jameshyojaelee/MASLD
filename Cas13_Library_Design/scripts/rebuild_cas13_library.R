#!/usr/bin/env Rscript
# rebuild_cas13_library.R
# ---------------------------------------------------------------------------
# Build cas13_library_v3.0.csv under the v6 definition (2026-06-04). v6 extends
# v5 (PI 2026-06-01, 6-team red-team) with: (1) the human spine sourced from
# limma-voom + metafor (ashr-shrunk) instead of dream; (2) a MASH-vs-MASL
# additive tier; (3) a miRNA tier. HUMAN-ANCHORED, READOUT-AWARE, REPRODUCIBLE.
#
#   LIBRARY = HUMAN SPINE  (mouse orthologs of human limma-voom+metafor ashr DEGs,
#                           lfsr<0.05 & shrunk_logFC>0.2, UP only)
#             UNION
#             MOUSE-CONFIRMED TIER (mouse cross-diet UP, ashr lfsr<0.05 &
#                           shrunk_logFC>0.5 in >=3 of 4 diets, AND the gene's
#                           human ortholog has human logFC>0 -- i.e. mouse
#                           evidence is admitted only with human directional
#                           concordance, never mouse-alone)
#             UNION
#             MASH-PROGRESSION TIER (mouse orthologs of MASH-vs-MASL UP human
#                           DEGs, lfsr<0.05 & shrunk_logFC>0.2; full additive)
#             UNION
#             miRNA TIER    (ortholog-conserved mouse<->human miRNAs, tier-H
#                           miRBase/MirGeneDB; DE-independent -- bulk polyA and
#                           10x scRNA both miss mature miRNAs)
#   restricted to protein_coding + lncRNA + miRNA biotypes.
#
# v6 NOTES:
#  - Human DEG canonical (LIBRARY SCOPE) is now limma-voom per-study + metafor
#    REML (06_meta_analysis.R) ashr-shrunk by 06b_meta_ashr_shrinkage.R. The
#    multi-evidence atlas / paper-wide canonical remains dream unless propagated.
#  - tier priority on overlap: human > mouse_confirmed > mash_progression >
#    mirna_conserved. is_mash_deg is an annotation on ALL rows (overlap-aware).
#  - miRNA hep-expression uses a miRNA-appropriate reference (mirna_hep_expression.csv:
#    liver small-RNA atlas + intragenic host-gene inheritance), since the scRNA
#    hep_substrate tag cannot score miRNAs. miRNAs are SOFT-tagged, not hard-dropped.
#  - Cas13 targets the pri-miRNA / host transcript for miRNA rows (mature ~22nt is
#    below the guide footprint) -- recorded as mirna_target_substrate="pri-miRNA".
#
# WHY v5 (changes from v4, all from the red-team):
#  - HUMAN-ANCHORED, not a flat union. The human arm is the unconditional spine;
#    mouse cross-diet genes enter only with human directional concordance. This
#    removes the 18% mouse-only tail that carried no human support.
#  - The hard 70%-human GATE is DELETED. In v4 the diet-replication threshold
#    (MIN_DIETS) was reverse-engineered to clear that gate (>=2 gave 65%<70%,
#    so >=3 was chosen). %human is now a REPORTED OUTCOME, not a constraint.
#  - Human effect-size floor 0.2 -> 0.3 (drops the weak 0.2-0.3 tail, ~1.23x,
#    where Cas13 knockdown has the least phenotypic headroom).
#  - SOFT readout-aware tags (NO filtering): hep_substrate (is the transcript a
#    valid hepatocyte Cas13 substrate, from per-lineage pseudobulk specificity)
#    and sc_disease_celltype (which cell type the disease-DE concentrates in).
#    The screen reads out CELL-AUTONOMOUS HEPATOCYTE lipid; ~60% of bulk-MASH
#    DEGs are non-hepatocyte. Tags make scoreability explicit; the pooled screen
#    tolerates non-scorers as internal negatives.
#  - REPRODUCIBILITY: deterministic ortholog tie-breaks (one2one-preferred, then
#    by ensembl id); n_human_orthologs / ortholog_ambiguous / direction_conflict
#    flags (fixes the v4 ACACA-row-labelled-ACACB direction-incoherence); a
#    provenance assert on the canonical STAR -s 2 ashr input; a BUILD_MANIFEST.
#
# NOTE on the earlier "stale -s 0 input" alarm: it was FALSE. The consumed
# dream_results_ashr.csv IS the canonical STAR -s 2 file (~27,638-gene universe);
# the misleadingly named *_star.csv is the OLD -s 0 run. The assert below guards it.
#
# COLOC/genetics stays an ANNOTATION flag (has_coloc, hep_expressed_causal), NOT
# an entry axis. Project "no discordance exclusion" principle preserved: the human
# SPINE admits human-up genes regardless of mouse direction; the human-concordance
# requirement applies ONLY to the mouse-confirmed secondary tier.
#
# Consumed by figS_cas13_library_*.R (S_lib_2..7). Output filename kept as
# cas13_library_v3.0.csv for consumer compatibility; library_version column = v5.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PERDIET <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
OUT     <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0.csv")
BACKUP  <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0_pre_v6.csv")
DIFFOUT <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v5_to_v6_diff.csv")
MANIFEST<- file.path(BASE, "Cas13_Library_Design/data/BUILD_MANIFEST.txt")
ATLAS   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
ORTHO   <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
META    <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
# v6: human DEG source is the limma-voom + metafor REML meta-analysis, ashr-shrunk
# (06b_meta_ashr_shrinkage.R), replacing dream. Schema matches dream_results_ashr.csv:
# gene, logFC, shrunk_logFC, lfsr, symbol. dream remains canonical for the atlas/paper.
ASHR    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/meta_results_ashr.csv")
# v6: MASH-vs-MASL (NASH>NAFL) limma-voom+metafor, ashr-shrunk (06b). Additive tier.
MASH    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nafl_vs_nash_meta_ashr.csv")
HEPSPEC <- file.path(BASE, "Cas13_Library_Design/data/hep_specificity.csv")
# v6: miRNA hepatic-expression annotation (liver small-RNA atlas + intragenic
# host-gene inheritance; scRNA hep_substrate is blind to miRNAs). Built externally.
MIRNAHEP<- file.path(BASE, "Cas13_Library_Design/data/mirna_hep_expression.csv")

DIETS            <- c("MCD", "CDAHFD", "Western", "HFD")
MIN_DIETS        <- 3L       # mouse cross-diet replication (default; sensitivity 2/3/4 reported)
LFSR_THR         <- 0.05     # mouse arm
SHRUNK_LFC_THR   <- 0.5      # mouse arm (stricter bar; mouse-only effect size)
HUMAN_LFSR_THR   <- 0.05     # human arm
HUMAN_SHRUNK_THR <- 0.2      # human + MASH arms (0.2 = F1-calibrated optimum; PI 2026-06-01)
COLOC_PP4_THR    <- 0.5      # annotation only
KEEP_BIOTYPES    <- c("protein_coding", "lncRNA", "miRNA")
LIB_VERSION      <- "v6"

strip_v <- function(x) sub("[.][0-9]+$", "", x)

# --- 0. Back up the existing (v4) library before overwriting -----------------
if (file.exists(OUT)) {
  file.copy(OUT, BACKUP, overwrite = TRUE)
  cat("backed up v4 ->", BACKUP, "\n")
}

# --- 1. Mouse cross-diet membership (ashr-shrunk effect size) ----------------
up_by_diet <- lapply(DIETS, function(d) {
  dt <- fread(file.path(PERDIET, paste0(d, "_de_results.csv")))
  dt[, gene_base := strip_v(gene)]
  dt[lfsr < LFSR_THR & shrunk_logFC > SHRUNK_LFC_THR, gene_base]
})
names(up_by_diet) <- DIETS

all_up <- sort(unique(unlist(up_by_diet)))
cross <- data.table(gene_id_mouse = all_up)
for (d in DIETS) cross[, (d) := gene_id_mouse %in% up_by_diet[[d]]]
cross[, n_diets_up := rowSums(.SD), .SDcols = DIETS]
cross[, diet_list := apply(.SD, 1, function(r) paste(DIETS[as.logical(r)], collapse = ",")),
      .SDcols = DIETS]

# --- 2. Human arm: limma-voom+metafor ashr DEGs + provenance assert ----------
ash <- fread(ASHR, select = c("gene", "logFC", "shrunk_logFC", "lfsr", "symbol"))
ash[, hb := strip_v(gene)]
universe_n <- nrow(ash)
col1a1_lfc <- ash[symbol == "COL1A1", logFC][1]
cat(sprintf("PROVENANCE metafor-ashr: universe=%d  COL1A1 logFC=%.4f\n", universe_n, col1a1_lfc))
# metafor universe = genes in >=3 cohorts (~26,011); guard against the old -s 0
# dream run (34,453) or a mis-pointed file.
if (universe_n > 30000)
  stop(sprintf("metafor-ashr universe=%d > 30000 (looks like the -s 0 dream run). Expected metafor (~26,011). Check ASHR path/file.", universe_n))
if (is.na(col1a1_lfc) || col1a1_lfc < 0.5)
  stop(sprintf("COL1A1 metafor logFC=%.3f unexpected (expected ~1.34, up in disease). Check ASHR file.", col1a1_lfc))

human_logfc <- ash[, .(hb, human_logFC = logFC)]      # for mouse-tier concordance
human_up <- ash[!is.na(lfsr) & lfsr < HUMAN_LFSR_THR & shrunk_logFC > HUMAN_SHRUNK_THR]
cat(sprintf("Human metafor-ashr UP (lfsr<%.2f & shrunk_logFC>%.2f): %d genes\n",
            HUMAN_LFSR_THR, HUMAN_SHRUNK_THR, nrow(human_up)))

# --- 3. Ortholog table (deterministic, one2one-preferred) --------------------
ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl", "human_ensembl", "human_symbol",
                          "confidence_tier", "is_one2one"))
ortho[, gene_id_mouse := strip_v(mouse_ensembl)]
ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, o2o := is_one2one %in% c(TRUE, "True", "TRUE", "true")]
ortho[, one2one_rank := ifelse(o2o, 0L, 1L)]          # one2one preferred

# count of distinct human orthologs per mouse gene (ambiguity flag)
nho <- ortho[, .(n_human_orthologs = uniqueN(human_ensembl)), by = gene_id_mouse]

# (a) best HUMAN per mouse -- annotation map (deterministic tie-break by ensembl)
setorder(ortho, gene_id_mouse, trank, one2one_rank, human_ensembl)
m2h <- ortho[, .(m2h_human_ensembl = human_ensembl[1],
                 m2h_human_symbol  = human_symbol[1]), by = gene_id_mouse]

# (b) best MOUSE per human -- backbone map (deterministic tie-break by ensembl)
ortho_byhuman <- copy(ortho)
setorder(ortho_byhuman, human_ensembl, trank, one2one_rank, mouse_ensembl)
ortho_byhuman <- unique(ortho_byhuman, by = "human_ensembl")

# --- 4a. HUMAN SPINE: one mouse ortholog per human ashr DEG ------------------
hu <- merge(human_up[, .(human_ensembl = hb, human_symbol = symbol,
                         human_logFC = logFC, human_shrunk = shrunk_logFC)],
            ortho_byhuman[, .(human_ensembl, gene_id_mouse)], by = "human_ensembl")
# multiple human genes -> same mouse: keep the entry human with strongest evidence
setorder(hu, gene_id_mouse, -human_shrunk, human_ensembl)
spine <- hu[, .(entry_human_ensembl = human_ensembl[1],
                gene_symbol_human   = human_symbol[1],
                entry_human_logFC   = human_logFC[1]), by = gene_id_mouse]
mouse_h <- spine$gene_id_mouse
cat("Human spine (mouse orthologs of human metafor-ashr DEGs, 1/human):", length(mouse_h), "\n")

# --- 4a'. MASH-PROGRESSION TIER: mouse orthologs of MASH-vs-MASL UP DEGs -------
# Full additive: every MASH-up human DEG (lfsr<0.05 & shrunk_logFC>thr; positive =
# higher in NASH/MASH vs NAFL/MASL) contributes its mouse ortholog. Overlap with
# the disease-vs-control spine is dedup'd by the union below and surfaced via the
# is_mash_deg annotation.
mash <- fread(MASH, select = c("gene", "logFC", "shrunk_logFC", "lfsr", "symbol"))
mash[, hb := strip_v(gene)]
mash_up <- mash[!is.na(lfsr) & lfsr < HUMAN_LFSR_THR & shrunk_logFC > HUMAN_SHRUNK_THR]
cat(sprintf("MASH-vs-MASL UP (lfsr<%.2f & shrunk_logFC>%.2f): %d human genes\n",
            HUMAN_LFSR_THR, HUMAN_SHRUNK_THR, nrow(mash_up)))
mh <- merge(mash_up[, .(human_ensembl = hb, mash_shrunk_logFC = shrunk_logFC, mash_lfsr = lfsr)],
            ortho_byhuman[, .(human_ensembl, gene_id_mouse)], by = "human_ensembl")
setorder(mh, gene_id_mouse, -mash_shrunk_logFC, human_ensembl)  # strongest MASH effect/mouse gene
mash_tier <- unique(mh, by = "gene_id_mouse")
mouse_mash <- mash_tier$gene_id_mouse
cat("MASH-progression tier (mouse orthologs of MASH-up human DEGs):", length(mouse_mash), "\n")

# --- 4b. MOUSE-CONFIRMED TIER: cross-diet UP + human directional concordance --
mouse_confirmed_at <- function(min_diets) {
  mc_genes <- cross[n_diets_up >= min_diets, gene_id_mouse]
  mc <- merge(data.table(gene_id_mouse = mc_genes), m2h, by = "gene_id_mouse", all.x = TRUE)
  mc <- merge(mc, human_logfc, by.x = "m2h_human_ensembl", by.y = "hb", all.x = TRUE)
  mc[!is.na(human_logFC) & human_logFC > 0, gene_id_mouse]
}
mouse_confirmed <- mouse_confirmed_at(MIN_DIETS)
cat(sprintf("Mouse-confirmed tier (>=%d diets + human concordance): %d genes\n",
            MIN_DIETS, length(mouse_confirmed)))

# --- 4c. miRNA TIER: ortholog-conserved mouse<->human miRNAs (DE-independent) --
# Neither bulk polyA nor 10x scRNA captures mature miRNAs, so the miRNA tier enters
# by ortholog CONSERVATION, not differential expression. SUPER-CONFIDENT only:
# require >=2 of {miRBase family-ID, MirGeneDB, biomaRt/Ensembl-Compara}. Sequence-
# homology (MMseqs2/BLAST/TOGA) is uninformative for ~22nt mature / ~70nt hairpin
# miRNAs and contributes 0 here. >=2-of-3 keeps every canonical hepatic miRNA
# (miR-122/148a/21/192/22/143/let-7/451a...) that a single-source filter would drop.
mir <- fread(cmd = paste0("zcat ", ORTHO),
             select = c("mouse_ensembl", "human_ensembl", "human_symbol", "mouse_biotype",
                        "tier_H_mirbase", "tier_H_mirgenedb", "tier_H_biomart",
                        "mirbase_family", "mirbase_arm"))
mir[, gene_id_mouse := strip_v(mouse_ensembl)]
tier1 <- function(x) as.integer(x %in% c(1, "1", TRUE))
mir[, mirna_n_evidence := tier1(tier_H_mirbase) + tier1(tier_H_mirgenedb) + tier1(tier_H_biomart)]
mir <- mir[mouse_biotype == "miRNA" & mirna_n_evidence >= 2L]   # super-confident: >=2 of 3 sources
setorder(mir, gene_id_mouse, -mirna_n_evidence, human_ensembl)  # 1 row per mouse miRNA gene
mir <- unique(mir, by = "gene_id_mouse")
mouse_mir <- mir$gene_id_mouse
cat("miRNA conserved tier (super-confident, >=2 of miRBase/MirGeneDB/biomaRt):", length(mouse_mir), "\n")

# --- 5. Assemble library: union of 4 tiers, PC+lncRNA+miRNA only --------------
lib_genes <- sort(Reduce(union, list(mouse_h, mouse_confirmed, mouse_mash, mouse_mir)))
mem <- data.table(gene_id_mouse = lib_genes)
# tier priority on overlap: human > mouse_confirmed > mash_progression > mirna_conserved
mem[, tier := fifelse(gene_id_mouse %in% mouse_h,         "human",
              fifelse(gene_id_mouse %in% mouse_confirmed, "mouse_confirmed",
              fifelse(gene_id_mouse %in% mouse_mash,      "mash_progression",
                                                          "mirna_conserved")))]
mem[, has_human_de := gene_id_mouse %in% mouse_h]     # disease-vs-control spine membership
mem[, is_mash_deg  := gene_id_mouse %in% mouse_mash]  # MASH-vs-MASL annotation (overlap-aware)
cat("Library union (pre biotype filter):", nrow(mem), "genes\n")

# diet membership
mem <- merge(mem, cross[, .(gene_id_mouse, n_diets_up, diet_list)],
             by = "gene_id_mouse", all.x = TRUE)
mem[is.na(n_diets_up), n_diets_up := 0L]
mem[is.na(diet_list), diet_list := ""]

# MASH-vs-MASL effect (annotation on all rows that are MASH-up DEGs)
mem <- merge(mem, mash_tier[, .(gene_id_mouse, mash_shrunk_logFC, mash_lfsr)],
             by = "gene_id_mouse", all.x = TRUE)

# miRNA family / arm / evidence count (miRNA tier rows)
mem <- merge(mem, mir[, .(gene_id_mouse, mirbase_family, mirbase_arm, mirna_n_evidence)],
             by = "gene_id_mouse", all.x = TRUE)

# entry-human (spine) + m2h (annotation); ship the direction-coherent label
mem <- merge(mem, spine[, .(gene_id_mouse, gene_symbol_human, entry_human_ensembl)],
             by = "gene_id_mouse", all.x = TRUE)
mem <- merge(mem, m2h, by = "gene_id_mouse", all.x = TRUE)
mem[is.na(gene_symbol_human), gene_symbol_human := m2h_human_symbol]   # mouse-confirmed label
# human ensembl used for downstream annotation joins (spine: entry; mouse: m2h best)
mem[, human_ensembl_annot := ifelse(tier == "human", entry_human_ensembl, m2h_human_ensembl)]

# ortholog ambiguity + direction-conflict (ACACA/ACACB-class) flags
mem <- merge(mem, nho, by = "gene_id_mouse", all.x = TRUE)
mem[is.na(n_human_orthologs), n_human_orthologs := 0L]
mem[, ortholog_ambiguous := n_human_orthologs > 1L]
mem <- merge(mem, ash[, .(hb, m2h_logFC = logFC)],
             by.x = "m2h_human_ensembl", by.y = "hb", all.x = TRUE)
mem[, direction_conflict := tier == "human" & !is.na(m2h_human_symbol) &
      m2h_human_symbol != gene_symbol_human & !is.na(m2h_logFC) & m2h_logFC < 0]

# biotype + mouse symbol
meta <- fread(META)
setnames(meta, c("mouse_ensembl_base", "mouse_symbol_gtf", "mouse_biotype"),
         c("gene_id_mouse", "gene_symbol_mouse", "biotype"), skip_absent = TRUE)
mem <- merge(mem, meta[, .(gene_id_mouse, gene_symbol_mouse, biotype)],
             by = "gene_id_mouse", all.x = TRUE)
mem[grepl("protein_coding", biotype), biotype := "protein_coding"]
mem[grepl("lncRNA|lincRNA", biotype), biotype := "lncRNA"]
n_before <- nrow(mem)
mem <- mem[biotype %in% KEEP_BIOTYPES]
cat("Biotype filter PC+lncRNA:", n_before, "->", nrow(mem), "genes\n")

# --- 6. Readout-aware soft tags (hepatocyte substrate + disease cell type) ----
if (file.exists(HEPSPEC)) {
  hep <- fread(HEPSPEC)
  mem <- merge(mem, hep[, .(gene_symbol_human = gene_symbol, hep_mean_cpm,
                            hep_ratio, hep_substrate, hep_top_other = top_other_celltype)],
               by = "gene_symbol_human", all.x = TRUE)
  mem[is.na(hep_substrate), hep_substrate := "unknown"]
} else {
  warning("hep_specificity.csv not found; run compute_hep_specificity.R first. Tagging 'unknown'.")
  mem[, `:=`(hep_mean_cpm = NA_real_, hep_ratio = NA_real_,
             hep_substrate = "unknown", hep_top_other = NA_character_)]
}

# disease-DE cell type + COLOC annotation from the multi-evidence atlas
atlas_sc <- fread(ATLAS, select = c("ensembl_id", "sc_best_celltype",
                                    "coloc_best_susie_pp4_polyfun"))
atlas_sc[, hb := strip_v(ensembl_id)]
mem <- merge(mem, atlas_sc[, .(hb, sc_disease_celltype = sc_best_celltype,
                               coloc_pp4 = coloc_best_susie_pp4_polyfun)],
             by.x = "human_ensembl_annot", by.y = "hb", all.x = TRUE)
mem[, has_coloc := !is.na(coloc_pp4) & coloc_pp4 > COLOC_PP4_THR]
mem[, hep_expressed_causal := has_coloc & hep_substrate == "high"]
mem[, has_human_evidence := has_human_de | has_coloc]

# lncRNA flags (orthology + hepatocyte expression)
mem[, lnc_human_ortholog := biotype == "lncRNA" & !is.na(gene_symbol_human) &
      gene_symbol_human != "" & !grepl("^ENSG", gene_symbol_human)]
mem[, lnc_hep_expressed := biotype == "lncRNA" & hep_substrate %in% c("high", "ambient_suspect")]

# miRNA hepatic-expression tags (miRNA-appropriate: liver small-RNA atlas +
# intragenic host-gene inheritance; the scRNA hep_substrate tag is blind to miRNAs).
# Soft-tag -- conserved miRNAs are kept regardless; the flag prioritizes scoreable ones.
if (file.exists(MIRNAHEP)) {
  mhep <- fread(MIRNAHEP)
  mhep[, gene_id_mouse := strip_v(mouse_ensembl)]
  mem <- merge(mem, unique(mhep[, .(gene_id_mouse, mirna_hep_expressed, mirna_hep_source)],
                           by = "gene_id_mouse"),
               by = "gene_id_mouse", all.x = TRUE)
} else {
  warning("mirna_hep_expression.csv not found; miRNA hep tags pending (run the miRNA atlas build, then re-run).")
  mem[, `:=`(mirna_hep_expressed = NA, mirna_hep_source = NA_character_)]
}
mem[biotype == "miRNA" & is.na(mirna_hep_source),    mirna_hep_source    := "none"]
mem[biotype == "miRNA" & is.na(mirna_hep_expressed), mirna_hep_expressed := FALSE]
# Cas13 targets the pri-miRNA / host transcript (mature ~22nt < guide footprint).
mem[, mirna_target_substrate := fifelse(biotype == "miRNA", "pri-miRNA", NA_character_)]

mem[, library_version := LIB_VERSION]

# --- 7. Write canonical schema (old columns preserved + v5 additions) --------
out <- mem[, .(gene_id_mouse, gene_symbol_mouse, gene_symbol_human, biotype, tier,
               n_diets_up, diet_list, has_human_de, is_mash_deg, mash_shrunk_logFC, mash_lfsr,
               has_coloc, has_human_evidence,
               hep_substrate, hep_mean_cpm, hep_ratio, sc_disease_celltype,
               n_human_orthologs, ortholog_ambiguous, direction_conflict,
               hep_expressed_causal, lnc_human_ortholog, lnc_hep_expressed,
               mirbase_family, mirbase_arm, mirna_n_evidence, mirna_hep_expressed, mirna_hep_source,
               mirna_target_substrate, library_version)]
setorder(out, tier, -n_diets_up, gene_symbol_mouse)
fwrite(out, OUT)

# --- 8. Verification + reported outcomes (NO hard gate) ----------------------
n_total <- nrow(out)
n_pc    <- out[biotype == "protein_coding", .N]
n_lnc   <- out[biotype == "lncRNA", .N]
n_mir   <- out[biotype == "miRNA", .N]
n_hd    <- out[tier == "human", .N]
pct_h   <- 100 * n_hd / n_total
n_sgrna <- (n_pc + n_lnc + n_mir) * 4L + 100L * 4L + 500L   # 4 gRNA/target (PC + lncRNA + miRNA)
# Coverage with delivery/sort efficiencies (PI 2026-06-01):
#   lenti transduction 30%; Cre recombination 80% (only Cre+ cells have active
#   RfxCas13d -> informative); FACS sorting efficiency 60% (sorted top/bottom
#   arms only; the unsorted input aliquot is set aside pre-sort, no FACS loss).
#   Effective informative cells/mouse = the smallest sequenced pool.
HEP <- 1e7; TRANSD <- 0.30; CRE_EFF <- 0.80; FACS_EFF <- 0.60
IN_FRAC <- 0.13; GATE <- 0.15; MORTALITY <- 1/8   # 1 in 8 mice die after viral injection
input_cells <- HEP * TRANSD * CRE_EFF * IN_FRAC                          # unsorted input (no FACS loss)
arm_cells   <- HEP * TRANSD * CRE_EFF * (1 - IN_FRAC) * GATE * FACS_EFF  # top/bottom 15% gate (FACS loss)
eff_cells   <- min(input_cells, arm_cells)
mice_cov <- ceiling(n_sgrna * 500 / eff_cells)             # SURVIVING mice needed for 500x
mice     <- ceiling(mice_cov / (1 - MORTALITY))            # INJECT this many to net the survivors

cat("\n=== LIBRARY", LIB_VERSION, "(", OUT, ") ===\n")
cat(sprintf("TOTAL: %d  (PC=%d, lncRNA=%d, miRNA=%d)\n", n_total, n_pc, n_lnc, n_mir))
cat("tier distribution (priority human>mouse_confirmed>mash_progression>mirna_conserved):\n")
print(table(out$tier))
cat(sprintf("  human spine = %.1f%% (REPORTED, not gated)\n", pct_h))
cat(sprintf("MASH-vs-MASL: is_mash_deg=%d (%.1f%% of library) | net-new mash_progression rows=%d\n",
            sum(out$is_mash_deg), 100 * mean(out$is_mash_deg), out[tier == "mash_progression", .N]))
if (n_mir > 0) {
  cat(sprintf("miRNA tier: %d genes | hep_expressed=%d\n",
              n_mir, out[biotype == "miRNA" & mirna_hep_expressed == TRUE, .N]))
  cat("  miRNA hep_source:\n"); print(table(out[biotype == "miRNA", mirna_hep_source]))
}
cat("hep_substrate distribution:\n"); print(table(out$hep_substrate))
cat(sprintf("scoreable core (hep_substrate==high): %d (%.1f%%)\n",
            out[hep_substrate == "high", .N], 100 * out[hep_substrate == "high", .N] / n_total))
cat("disease-DE cell type (top 6):\n"); print(head(sort(table(out$sc_disease_celltype), decreasing = TRUE), 6))
cat(sprintf("lncRNA: %d total | human-orthologous=%d | hep-expressed=%d\n",
            n_lnc, out[lnc_human_ortholog == TRUE, .N], out[lnc_hep_expressed == TRUE, .N]))
cat(sprintf("ortholog_ambiguous=%d | direction_conflict=%d | hep_expressed_causal=%d | has_coloc=%d\n",
            sum(out$ortholog_ambiguous), sum(out$direction_conflict),
            sum(out$hep_expressed_causal), sum(out$has_coloc)))
cat(sprintf("sgRNAs: %d (4/PC, 4/lncRNA, +100 ctrl x4, +500 NT)\n", n_sgrna))
cat(sprintf("coverage: transduction=%.0f%% Cre=%.0f%% FACS=%.0f%% | binding pool=%.0f cells/mouse (%s)\n",
            TRANSD*100, CRE_EFF*100, FACS_EFF*100, eff_cells,
            ifelse(arm_cells <= input_cells, "sorted arm", "input")))
cat(sprintf("mice: %d surviving needed @500x; INJECT %d (1/8=%.1f%% post-injection mortality)\n",
            mice_cov, mice, MORTALITY*100))

# ACACA spot-check (the v4 direction-incoherence)
acaca <- out[gene_symbol_mouse == "Acaca"]
if (nrow(acaca)) {
  cat(sprintf("ACACA check: gene_symbol_human=%s direction_conflict=%s\n",
              acaca$gene_symbol_human[1], acaca$direction_conflict[1]))
}

# MIN_DIETS sensitivity (human spine constant; mouse-confirmed varies)
cat("\nMIN_DIETS sensitivity (library size | %human):\n")
for (md in c(2L, 3L, 4L)) {
  mc <- mouse_confirmed_at(md)
  g  <- union(mouse_h, mc)
  bt <- meta[gene_id_mouse %in% g & grepl("protein_coding|lncRNA|lincRNA", biotype), gene_id_mouse]
  n  <- length(bt)
  nh <- length(intersect(bt, mouse_h))
  cat(sprintf("  >=%d diets: %d genes | %.1f%% human\n", md, n, 100 * nh / max(n, 1)))
}

# --- 9. v5 -> v6 diff --------------------------------------------------------
# Diff against the STABLE frozen v5 reference (cas13_library_v5_canonical.csv), so
# re-runs always show the true v5->v6 delta -- not the rolling BACKUP, which a
# re-run would overwrite with the previous v6 (giving a spurious empty diff).
V5REF   <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v5_canonical.csv")
diffref <- if (file.exists(V5REF)) V5REF else BACKUP
if (file.exists(diffref)) {
  v5 <- fread(diffref)
  v5g <- v5$gene_id_mouse; v6g <- out$gene_id_mouse
  added   <- setdiff(v6g, v5g)
  dropped <- setdiff(v5g, v6g)
  diff <- rbind(
    data.table(gene_id_mouse = added,   change = "added"),
    data.table(gene_id_mouse = dropped, change = "dropped")
  )
  diff <- merge(diff, out[, .(gene_id_mouse, gene_symbol_mouse, gene_symbol_human, tier, biotype)],
                by = "gene_id_mouse", all.x = TRUE)
  fwrite(diff, DIFFOUT)
  cat(sprintf("\nv5->v6 diff: +%d added / -%d dropped (v5=%d, v6=%d) -> %s\n",
              length(added), length(dropped), length(v5g), length(v6g), DIFFOUT))
}

# --- 10. BUILD_MANIFEST (provenance stamp) ----------------------------------
md5 <- function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
git_sha <- tryCatch(system2("git", c("-C", BASE, "rev-parse", "--short", "HEAD"),
                            stdout = TRUE, stderr = FALSE), error = function(e) "NA")
manifest <- c(
  sprintf("library_version: %s", LIB_VERSION),
  sprintf("built_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("git_sha: %s", paste(git_sha, collapse = "")),
  "human_deg_source: limma-voom per-study + metafor REML, ashr-shrunk (06b); dream retired for library",
  sprintf("n_genes: %d (PC=%d lncRNA=%d miRNA=%d)", n_total, n_pc, n_lnc, n_mir),
  sprintf("tier: human=%d mouse_confirmed=%d mash_progression=%d mirna_conserved=%d",
          out[tier == "human", .N], out[tier == "mouse_confirmed", .N],
          out[tier == "mash_progression", .N], out[tier == "mirna_conserved", .N]),
  sprintf("is_mash_deg (annotation, overlap-aware): %d", sum(out$is_mash_deg)),
  sprintf("pct_human_reported: %.1f", pct_h),
  sprintf("n_sgrna: %d  mice_surviving_500x: %d  mice_to_inject: %d", n_sgrna, mice_cov, mice),
  sprintf("coverage: transduction=%.2f Cre=%.2f FACS=%.2f IN_FRAC=%.2f GATE=%.2f mortality=%.3f -> eff_cells/mouse=%.0f (%s binds)",
          TRANSD, CRE_EFF, FACS_EFF, IN_FRAC, GATE, MORTALITY, eff_cells,
          ifelse(arm_cells <= input_cells, "arm", "input")),
  sprintf("HUMAN_SHRUNK_THR: %.2f  SHRUNK_LFC_THR: %.2f  MIN_DIETS: %d",
          HUMAN_SHRUNK_THR, SHRUNK_LFC_THR, MIN_DIETS),
  "inputs (path | mtime | md5):",
  sprintf("  human_metafor_ashr: %s | %s | %s", ASHR, format(file.mtime(ASHR)), md5(ASHR)),
  sprintf("  mash_metafor_ashr: %s | %s | %s", MASH, format(file.mtime(MASH)), md5(MASH)),
  sprintf("  ortho: %s | %s | %s", ORTHO, format(file.mtime(ORTHO)), md5(ORTHO)),
  sprintf("  atlas: %s | %s | %s", ATLAS, format(file.mtime(ATLAS)), md5(ATLAS)),
  sprintf("  hepspec: %s | %s | %s", HEPSPEC, format(file.mtime(HEPSPEC)), md5(HEPSPEC)),
  sprintf("  mirna_hep: %s | %s | %s", MIRNAHEP, format(file.mtime(MIRNAHEP)), md5(MIRNAHEP)),
  sprintf("  metafor_universe: %d  COL1A1_logFC: %.4f (metafor REML signature, up in disease)", universe_n, col1a1_lfc),
  paste0("  per_diet: ", paste(DIETS, collapse = ","))
)
writeLines(manifest, MANIFEST)
cat("wrote", MANIFEST, "\n")
