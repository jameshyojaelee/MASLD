#!/usr/bin/env Rscript
# Phase 0 hit-set builder for perturbation campaign (5 arms: D1-D5)
# Outputs to Analysis/Perturbation/data/hits/

suppressPackageStartupMessages({
  library(data.table)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HITS_DIR <- file.path(PROJ, "Analysis/Perturbation/data/hits")
dir.create(HITS_DIR, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                 paste0(...)))

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
msg("Loading convergence_evidence.csv")
ce <- fread(file.path(PROJ,
                      "RNA-seq/results/multi_evidence/convergence_evidence.csv"))
setnames(ce, "human_symbol", "gene")
msg("  ", nrow(ce), " rows")

msg("Loading multi_evidence_atlas.csv")
atlas <- fread(file.path(PROJ,
                         "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
setnames(atlas, "human_symbol", "gene")
msg("  ", nrow(atlas), " rows x ", ncol(atlas), " cols")

msg("Loading gene_level_coloc.csv")
coloc <- fread(file.path(PROJ,
                          "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
msg("  ", nrow(coloc), " rows")

msg("Loading convergent_drug_targets.csv")
drugs <- fread(file.path(PROJ,
                          "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv"))
msg("  ", nrow(drugs), " rows")

msg("Loading hepatocyte subtype_markers.csv + meta_subtype_mapping.csv")
hep_mk <- fread(file.path(PROJ,
            "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv"))
hep_meta <- fread(file.path(PROJ,
            "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv"))
hep_mk <- merge(hep_mk, hep_meta[, .(subtype, meta_subtype)],
                by = "subtype", all.x = TRUE)
msg("  ", nrow(hep_mk), " marker rows; meta classes: ",
    paste(unique(hep_mk$meta_subtype), collapse = ", "))

msg("Loading Conserved gene list (gene_concordance_per_gene.csv)")
gcp <- fread(file.path(PROJ,
            "Analysis/Cross_Species_Concordance/results/gene_concordance_per_gene.csv"))
conserved_core <- unique(gcp[primary_category == "Conserved",
                              .(mouse_gene_id, human_symbol)])
msg("  Conserved_Core: ", nrow(conserved_core), " genes")

msg("Loading mouse_human_orthologs.tsv.gz")
ortho <- fread(file.path(PROJ,
            "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"))
setnames(ortho,
         c("mouse_ensembl_gene_id", "human_ensembl_gene_id", "orthology_type"),
         c("mouse_ensembl", "human_ensembl", "orthology_type"))
msg("  ", nrow(ortho), " ortholog pairs")

# Build human_ensembl -> human_symbol map from atlas
gene_map <- atlas[!is.na(gene) & gene != "",
                  .(human_ensembl = ensembl_id, human_symbol = gene)]
gene_map <- unique(gene_map, by = "human_ensembl")

# Build mouse_ensembl -> mouse_symbol from conserved_core (only where available)
mouse_map <- unique(gcp[, .(mouse_gene_id, human_symbol)])

# -----------------------------------------------------------------------------
# D1 -- mechanism hit-set (~800 genes)
# -----------------------------------------------------------------------------
msg("Building D1 mechanism hit-set")
# (a) Top 500 from convergence (Tier 1 + Tier 2)
ce_t12 <- ce[!excluded_from_ranking & tier %in% c("Tier1", "Tier2", 1, 2, "1", "2") &
              !is.na(convergence_rank)]
if (nrow(ce_t12) == 0) {
  # tiers may be coded differently; fall back to top by rank
  ce_t12 <- ce[!excluded_from_ranking & !is.na(convergence_rank)]
}
ce_top <- ce_t12[order(convergence_rank)][1:min(500, .N)]
ce_top[, source_d1_convergence := TRUE]

# (b) COLOC PP4>0.8 (SuSiE PP4 preferred; fall back to ABF PP4)
coloc_d1 <- coloc[(!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > 0.8) |
                  (!is.na(coloc_best_pp4) & coloc_best_pp4 > 0.8)]
msg("  COLOC PP4>0.8: ", nrow(coloc_d1), " genes")

# (c) Top 200 drug targets by causal_score / convergence_score
drug_score_col <- intersect(c("causal_score", "convergence_score"),
                             names(drugs))[1]
if (is.na(drug_score_col) || is.null(drug_score_col)) drug_score_col <- "convergence_score"
drugs_d1 <- drugs[order(-get(drug_score_col))][1:min(200, .N)]
msg("  Top drug targets (by ", drug_score_col, "): ", nrow(drugs_d1))

# (d) Hepatocyte subtype markers: Progressor + Healthy (top by score within each subtype)
hep_progr <- hep_mk[meta_subtype == "Disease-Progressor"]
hep_health <- hep_mk[meta_subtype == "Healthy"]
# top 50 each (then dedup)
hep_progr_top <- unique(hep_progr[order(-scores), .(gene = names)])[1:min(50, .N)]
hep_health_top <- unique(hep_health[order(-scores), .(gene = names)])[1:min(50, .N)]
hep_top <- unique(rbind(hep_progr_top, hep_health_top))
msg("  Hepatocyte subtype markers (Progressor+Healthy top): ", nrow(hep_top))

# Assemble D1
d1 <- data.table(gene = character(), sources = character())

# Convergence
d1_add <- ce_top[, .(gene, source = "convergence_top500",
                     convergence_rank = as.numeric(convergence_rank),
                     coloc_pp4 = NA_real_, druggable = FALSE,
                     hep_subtype_marker = "")]
# Add COLOC
coloc_pp4_best <- pmax(coloc$coloc_best_susie_pp4, coloc$coloc_best_pp4,
                       na.rm = TRUE)
coloc_lookup <- coloc[, .(gene, coloc_pp4 = pmax(coloc_best_susie_pp4,
                                                  coloc_best_pp4, na.rm = TRUE))]
coloc_lookup <- coloc_lookup[!is.na(coloc_pp4)]
setorder(coloc_lookup, -coloc_pp4)
coloc_lookup <- unique(coloc_lookup, by = "gene")

d1_coloc <- merge(coloc_d1[, .(gene)],
                  coloc_lookup, by = "gene", all.x = TRUE)
d1_coloc[, `:=`(source = "coloc_pp4_gt_0.8",
                convergence_rank = NA_real_,
                druggable = FALSE,
                hep_subtype_marker = "")]

# Drug
d1_drugs <- drugs_d1[, .(gene = symbol, source = "drug_target",
                          convergence_rank = NA_real_,
                          coloc_pp4 = NA_real_, druggable = TRUE,
                          hep_subtype_marker = "")]

# Hep
d1_hep_progr_dt <- hep_progr_top[, .(gene, source = "hep_marker",
                                     convergence_rank = NA_real_,
                                     coloc_pp4 = NA_real_, druggable = FALSE,
                                     hep_subtype_marker = "Disease-Progressor")]
d1_hep_health_dt <- hep_health_top[, .(gene, source = "hep_marker",
                                       convergence_rank = NA_real_,
                                       coloc_pp4 = NA_real_, druggable = FALSE,
                                       hep_subtype_marker = "Healthy")]

d1_all <- rbind(d1_add[, .(gene, source, convergence_rank, coloc_pp4, druggable,
                            hep_subtype_marker)],
                d1_coloc[, .(gene, source, convergence_rank, coloc_pp4, druggable,
                              hep_subtype_marker)],
                d1_drugs, d1_hep_progr_dt, d1_hep_health_dt)
d1_all <- d1_all[!is.na(gene) & gene != ""]

# Aggregate sources per gene
d1_agg <- d1_all[, .(
  sources = paste(sort(unique(source)), collapse = ";"),
  convergence_rank = suppressWarnings(min(convergence_rank, na.rm = TRUE)),
  coloc_pp4 = suppressWarnings(max(coloc_pp4, na.rm = TRUE)),
  druggable = any(druggable),
  hep_subtype_marker = paste(sort(unique(hep_subtype_marker[hep_subtype_marker != ""])),
                              collapse = ";")
), by = gene]
d1_agg[is.infinite(convergence_rank), convergence_rank := NA_real_]
d1_agg[is.infinite(coloc_pp4), coloc_pp4 := NA_real_]
setorder(d1_agg, convergence_rank, -coloc_pp4)

fwrite(d1_agg, file.path(HITS_DIR, "d1_mechanism_hits.csv"))
msg("  D1 wrote ", nrow(d1_agg), " unique genes -> d1_mechanism_hits.csv")
n_d1 <- nrow(d1_agg)

# -----------------------------------------------------------------------------
# D2 -- reversal hit-set (all atlas genes)
# -----------------------------------------------------------------------------
msg("Building D2 reversal hit-set (all atlas genes)")
# Merge convergence_rank + coloc into atlas
atlas_d2 <- atlas[, .(gene, ensembl_id, gene_biotype)]
ce_lookup <- ce[, .(gene, convergence_rank)]
ce_lookup <- ce_lookup[!duplicated(gene)]
atlas_d2 <- merge(atlas_d2, ce_lookup, by = "gene", all.x = TRUE)
atlas_d2 <- merge(atlas_d2, coloc_lookup, by = "gene", all.x = TRUE)
# dream_logFC/padj already in atlas
dream_lookup <- atlas[, .(gene, dream_logFC, dream_padj)]
dream_lookup <- dream_lookup[!duplicated(gene)]
atlas_d2 <- merge(atlas_d2, dream_lookup, by = "gene", all.x = TRUE)
atlas_d2[, atlas_in := TRUE]
setcolorder(atlas_d2, c("gene", "ensembl_id", "gene_biotype", "atlas_in",
                         "convergence_rank", "coloc_pp4", "dream_logFC",
                         "dream_padj"))
fwrite(atlas_d2, file.path(HITS_DIR, "d2_reversal_hits.csv"))
msg("  D2 wrote ", nrow(atlas_d2), " genes -> d2_reversal_hits.csv")
n_d2 <- nrow(atlas_d2)

# -----------------------------------------------------------------------------
# D3 -- synergy hit-sets (multi-tier)
# -----------------------------------------------------------------------------
msg("Building D3 synergy hit-sets")

# Tier 1: top 100 convergence -> 100C2 = 4,950 pairs
top100 <- ce[!excluded_from_ranking & !is.na(convergence_rank)
             ][order(convergence_rank)][1:100]
top100_genes <- top100$gene
rank_lookup <- setNames(top100$convergence_rank, top100$gene)

pairs_t1 <- as.data.table(t(combn(top100_genes, 2)))
setnames(pairs_t1, c("V1", "V2"), c("gene1", "gene2"))
pairs_t1[, pair_source_tier := "tier1_top100_convergence"]
pairs_t1[, sum_convergence_rank := rank_lookup[gene1] + rank_lookup[gene2]]
fwrite(pairs_t1, file.path(HITS_DIR, "d3_tier1_pairs.csv"))
msg("  D3 tier1 pairs: ", nrow(pairs_t1), " -> d3_tier1_pairs.csv")
n_d3_t1 <- nrow(pairs_t1)

# Tier 2: all COLOC PP4>0.5 -> 662C2 = 218,791 pairs
coloc_pp5 <- coloc[(!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > 0.5) |
                   (!is.na(coloc_best_pp4) & coloc_best_pp4 > 0.5),
                   .(gene)]
coloc_pp5 <- unique(coloc_pp5[gene %in% atlas$gene])
msg("  COLOC PP4>0.5 unique genes (in atlas): ", nrow(coloc_pp5))
# Get convergence rank for ranking (NA -> Inf)
coloc_rank <- ce_lookup[gene %in% coloc_pp5$gene]
all_coloc_genes <- sort(coloc_pp5$gene)
big_rank_lookup <- setNames(rep(NA_integer_, length(all_coloc_genes)),
                             all_coloc_genes)
big_rank_lookup[coloc_rank$gene] <- coloc_rank$convergence_rank

# Generate pairs efficiently (combn is fine for 662)
pairs_mat <- combn(all_coloc_genes, 2)
pairs_t2 <- data.table(gene1 = pairs_mat[1, ], gene2 = pairs_mat[2, ])
pairs_t2[, pair_source_tier := "tier2_coloc_pp4_gt_0.5"]
pairs_t2[, sum_convergence_rank := big_rank_lookup[gene1] +
            big_rank_lookup[gene2]]
fwrite(pairs_t2, file.path(HITS_DIR, "d3_tier2_pairs.csv"))
msg("  D3 tier2 pairs: ", nrow(pairs_t2), " -> d3_tier2_pairs.csv")
n_d3_t2 <- nrow(pairs_t2)

# Tier 3: top 50 -> 50C3 = 19,600 triples
top50 <- ce[!excluded_from_ranking & !is.na(convergence_rank)
            ][order(convergence_rank)][1:50]
top50_rank <- setNames(top50$convergence_rank, top50$gene)
triples_mat <- combn(top50$gene, 3)
triples_t3 <- data.table(gene1 = triples_mat[1, ],
                          gene2 = triples_mat[2, ],
                          gene3 = triples_mat[3, ])
triples_t3[, pair_source_tier := "tier3_top50_convergence"]
triples_t3[, sum_convergence_rank := top50_rank[gene1] + top50_rank[gene2] +
             top50_rank[gene3]]
fwrite(triples_t3, file.path(HITS_DIR, "d3_tier3_triples.csv"))
msg("  D3 tier3 triples: ", nrow(triples_t3), " -> d3_tier3_triples.csv")
n_d3_t3 <- nrow(triples_t3)

# Tier 4: top 20 -> 20C4 = 4,845 quadruples
top20 <- ce[!excluded_from_ranking & !is.na(convergence_rank)
            ][order(convergence_rank)][1:20]
top20_rank <- setNames(top20$convergence_rank, top20$gene)
quads_mat <- combn(top20$gene, 4)
quads_t4 <- data.table(gene1 = quads_mat[1, ], gene2 = quads_mat[2, ],
                        gene3 = quads_mat[3, ], gene4 = quads_mat[4, ])
quads_t4[, pair_source_tier := "tier4_top20_convergence"]
quads_t4[, sum_convergence_rank := top20_rank[gene1] + top20_rank[gene2] +
            top20_rank[gene3] + top20_rank[gene4]]
fwrite(quads_t4, file.path(HITS_DIR, "d3_tier4_quadruples.csv"))
msg("  D3 tier4 quadruples: ", nrow(quads_t4), " -> d3_tier4_quadruples.csv")
n_d3_t4 <- nrow(quads_t4)

# Tier 5 placeholder
fwrite(data.table(gene1 = character(), gene2 = character(),
                   pair_source_tier = character(),
                   sum_convergence_rank = numeric(),
                   note = character()),
       file.path(HITS_DIR, "d3_tier5_adaptive.csv"))
msg("  D3 tier5 placeholder written -> d3_tier5_adaptive.csv")

# -----------------------------------------------------------------------------
# D4 -- circuit hit-set (hepatocyte ligands across receivers)
# -----------------------------------------------------------------------------
msg("Building D4 circuit hit-set")
d4_files <- c(
  Macrophages       = file.path(PROJ,
        "Analysis/SingleCell/results_gpu_v2/ccc/D2_hep_Mac_bidirectional.csv"),
  Stellate          = file.path(PROJ,
        "Analysis/SingleCell/results_gpu_v2/ccc/D2_hep_HSC_bidirectional.csv"),
  LSEC              = file.path(PROJ,
        "Analysis/SingleCell/results_gpu_v2/ccc/D2_hep_LSEC_bidirectional.csv")
)
# Note: Cholangiocyte file (D2_hep_Chol_*) was not produced upstream;
# only 3 of 4 receivers are available. Documented in README.
d4_list <- list()
for (rec in names(d4_files)) {
  f <- d4_files[[rec]]
  if (!file.exists(f)) {
    msg("  WARN: ", f, " missing")
    next
  }
  d <- fread(f)
  # Keep Hep_to_<receiver> direction only -> hepatocyte ligand
  d_lig <- d[grepl("^Hep_to_", direction)]
  setorder(d_lig, -score_diff)
  d_lig <- d_lig[1:min(50, .N)]
  d_lig[, receiver_cell_type := rec]
  d_lig[, source_file := basename(f)]
  d4_list[[rec]] <- d_lig[, .(ligand, receiver_cell_type, receptor,
                                ccc_score = score_diff, source_file)]
}
d4 <- rbindlist(d4_list)
fwrite(d4, file.path(HITS_DIR, "d4_circuit_ligands.csv"))
msg("  D4 wrote ", nrow(d4), " ligand-receiver-receptor rows -> d4_circuit_ligands.csv")
n_d4 <- nrow(d4)

# -----------------------------------------------------------------------------
# D5 -- mouse ortholog hit-set
# -----------------------------------------------------------------------------
msg("Building D5 mouse ortholog hit-set")
# Conserved_Core: 1,108 human-symbol genes (already have mouse_gene_id)
d5_conserved <- unique(gcp[primary_category == "Conserved",
                            .(human_symbol, mouse_ensembl = mouse_gene_id,
                              conserved_core = TRUE)])

# Top 100 human convergence
top100_d5 <- ce[!excluded_from_ranking & !is.na(convergence_rank)
                ][order(convergence_rank)][1:100,
                                            .(human_symbol = gene,
                                              convergence_rank)]
d5_conv <- merge(top100_d5, unique(gene_map[, .(human_symbol, human_ensembl)]),
                  by = "human_symbol", all.x = TRUE)
# Map via ortholog file (strict 1:1 vs other)
d5_conv <- merge(d5_conv, ortho, by = "human_ensembl", all.x = TRUE)
d5_conv[, conserved_core := FALSE]

# Conserved core: try to map mouse_ensembl -> ortholog_type from ortho
d5_conserved2 <- merge(d5_conserved, ortho[, .(mouse_ensembl, orthology_type)],
                        by = "mouse_ensembl", all.x = TRUE)
d5_conserved2 <- merge(d5_conserved2,
                        ce_lookup[, .(human_symbol = gene, convergence_rank)],
                        by = "human_symbol", all.x = TRUE)

# Merge into common schema
d5_a <- d5_conserved2[, .(human_gene = human_symbol,
                           mouse_ensembl,
                           ortholog_type = orthology_type,
                           conserved_core = TRUE,
                           convergence_rank)]
d5_b <- d5_conv[, .(human_gene = human_symbol,
                     mouse_ensembl,
                     ortholog_type = orthology_type,
                     conserved_core = FALSE,
                     convergence_rank)]
# Drop duplicates: if a top-100 conv gene is already in Conserved_Core, keep conserved flag
d5_b <- d5_b[!human_gene %in% d5_a$human_gene]
d5_full <- rbind(d5_a, d5_b)

# Attach mouse symbol from gcp (mouse_gene_id, human_symbol -> mouse_symbol via celltype file)
ccg <- fread(file.path(PROJ,
              "Analysis/Cross_Species_Concordance/results/celltype_conserved_genes_all.csv"))
mouse_sym_map <- unique(ccg[, .(mouse_ensembl = mouse_gene, mouse_gene = mouse_symbol)])
d5_full <- merge(d5_full, mouse_sym_map, by = "mouse_ensembl", all.x = TRUE)

# Flag missing
d5_full[is.na(ortholog_type) | ortholog_type == "",
         ortholog_type := ifelse(is.na(mouse_ensembl), "missing", "unknown")]
d5_full[is.na(mouse_ensembl), mouse_ensembl := ""]
d5_full[is.na(mouse_gene), mouse_gene := ""]

setcolorder(d5_full, c("human_gene", "mouse_gene", "ortholog_type",
                         "conserved_core", "convergence_rank", "mouse_ensembl"))
setorder(d5_full, -conserved_core, convergence_rank)
fwrite(d5_full, file.path(HITS_DIR, "d5_mouse_orthologs.csv"))
msg("  D5 wrote ", nrow(d5_full), " human-mouse pairs -> d5_mouse_orthologs.csv")
n_d5 <- nrow(d5_full)
n_d5_one2one <- sum(d5_full$ortholog_type == "ortholog_one2one")
n_d5_missing <- sum(d5_full$ortholog_type %in% c("missing", "", "unknown"))

# -----------------------------------------------------------------------------
# README
# -----------------------------------------------------------------------------
readme <- sprintf("# Perturbation Campaign -- Phase 0 hit-sets

Built %s by `scripts/build_hit_sets.R` (phase0-hitset-builder).
Per-arm \"what to perturb\" input lists for the Model Runners.

Source inputs (read-only):
- `RNA-seq/results/multi_evidence/multi_evidence_atlas.csv` (33,943 genes x 240 cols, 2026-05-19)
- `RNA-seq/results/multi_evidence/convergence_evidence.csv` (Tier + concordance + rank)
- `GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv` (canonical SuSiE+ABF gene-level)
- `RNA-seq/results/drug_repurposing/convergent_drug_targets.csv` (284 targets)
- `Analysis/SingleCell/results_gpu_v2/ccc/D2_hep_{Mac,HSC,LSEC}_bidirectional.csv`
- `Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv` + `meta_subtype_mapping.csv`
- `Analysis/Cross_Species_Concordance/results/gene_concordance_per_gene.csv` (Conserved_Core list)
- `archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz`

## Files (rows)

| File | Rows | Description |
|------|------|-------------|
| `d1_mechanism_hits.csv` | %d | D1 mechanism (top500 convergence + COLOC PP4>0.8 + top200 drug targets + Progressor/Healthy hep markers) |
| `d2_reversal_hits.csv` | %d | D2 reversal: all atlas genes (genome-wide) with biotype + convergence_rank + COLOC + dream LFC/padj |
| `d3_tier1_pairs.csv` | %d | D3 tier 1: top 100 convergence -> 100C2 pairs |
| `d3_tier2_pairs.csv` | %d | D3 tier 2: COLOC PP4>0.5 -> nC2 pairs |
| `d3_tier3_triples.csv` | %d | D3 tier 3: top 50 convergence -> 50C3 triples |
| `d3_tier4_quadruples.csv` | %d | D3 tier 4: top 20 convergence -> 20C4 quadruples |
| `d3_tier5_adaptive.csv` | 0 | D3 tier 5 placeholder; populated later from tier 1/2 screen results |
| `d4_circuit_ligands.csv` | %d | D4 hepatocyte ligands x 3 receivers (Macrophages, Stellate, LSEC) |
| `d5_mouse_orthologs.csv` | %d | D5 mouse ortholog map (Conserved_Core + top 100 convergence); %d one2one, %d missing |

## D1 mechanism schema

`gene`, `sources` (semicolon-sep), `convergence_rank`, `coloc_pp4`, `druggable`, `hep_subtype_marker`.
sources value vocabulary: `convergence_top500`, `coloc_pp4_gt_0.8`, `drug_target`, `hep_marker`.

## D2 reversal schema

`gene`, `ensembl_id`, `gene_biotype`, `atlas_in`, `convergence_rank`, `coloc_pp4`, `dream_logFC`, `dream_padj`.
Includes all biotypes (protein-coding + lncRNA + pseudogene + other); downstream filters with `gene_biotype`.

## D3 synergy schema

`gene1`, `gene2`, (`gene3`, `gene4` for tier3/4), `pair_source_tier`, `sum_convergence_rank`.
Tier-2 csv is the largest (~218k rows); consider streaming reader downstream.

## D4 circuit schema

`ligand`, `receiver_cell_type` (Macrophages / Stellate / LSEC), `receptor`, `ccc_score` (= score_diff
MASLD - control from D2_hep_*_bidirectional.csv `Hep_to_*` rows), `source_file`.

**Cholangiocyte gap**: upstream `Analysis/SingleCell/results_gpu_v2/ccc/` produced 3 of the 4 planned
hepatocyte->receiver paracrine files (Mac, HSC, LSEC) but no cholangiocyte file. The plan's 4-receiver
spec is reduced to 3 here; cholangiocyte CCC requires a new run by the CCC pipeline owner.

## D5 mouse ortholog schema

`human_gene`, `mouse_gene`, `ortholog_type` (ortholog_one2one / ortholog_one2many / ortholog_many2many /
missing / unknown), `conserved_core` (TRUE for the 1,108-gene Conserved_Core), `convergence_rank`,
`mouse_ensembl`.

## Re-run

```bash
micromamba activate rnaseq
Rscript Analysis/Perturbation/scripts/build_hit_sets.R
```
",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  n_d1, n_d2, n_d3_t1, n_d3_t2, n_d3_t3, n_d3_t4, n_d4, n_d5,
  n_d5_one2one, n_d5_missing
)

writeLines(readme, file.path(HITS_DIR, "README.md"))
msg("Wrote README.md")

msg("DONE -- ", HITS_DIR)
