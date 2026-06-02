#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_268_enrichment
#SBATCH --output=logs/net_268_enrichment_%j.out
#SBATCH --error=logs/net_268_enrichment_%j.err
# ===========================================================================
# Script 268: Community enrichment and biological labeling
# ===========================================================================
# Purpose: For each macro and meso community, run ORA (Fisher's exact test)
#          against MSigDB Hallmark pathways, drug targets, COLOC genes, sex
#          class, progression genes, and zonation class.  Assign a biological
#          label to each community based on the top enriched pathway.
#
# Input:
#   - RNA-seq/results/network/communities/community_assignments.csv
#   - RNA-seq/results/network/network_nodes.csv
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Output:
#   - RNA-seq/results/network/communities/community_enrichment.csv
#   - RNA-seq/results/network/communities/community_labels.csv
#   - RNA-seq/results/network/communities/communities.json
#
# Environment: rnaseq (data.table, msigdbr v10, jsonlite)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(jsonlite)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NET_DIR  <- file.path(BASE, "RNA-seq/results/network")
COMM_DIR <- file.path(NET_DIR, "communities")
ATLAS_PATH <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

cat("=== 268: Community Enrichment and Labeling ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ======================================================================
# 1. Load data
# ======================================================================
cat("--- Loading data ---\n")

# Community assignments
comm_path <- file.path(COMM_DIR, "community_assignments.csv")
if (!file.exists(comm_path)) {
  stop("Community assignments not found: ", comm_path, "\n  Run Script 266 first.")
}
comm <- fread(comm_path)
cat(sprintf("  Communities: %d genes\n", nrow(comm)))

# Node set
node_path <- file.path(NET_DIR, "network_nodes.csv")
if (!file.exists(node_path)) {
  stop("Node file not found: ", node_path)
}
nodes <- fread(node_path)
gene_col <- ifelse("gene" %in% names(nodes), "gene", "human_symbol")
cat(sprintf("  Nodes: %d genes\n", nrow(nodes)))

# Atlas
if (!file.exists(ATLAS_PATH)) {
  stop("Atlas not found: ", ATLAS_PATH)
}
atlas <- fread(ATLAS_PATH)
cat(sprintf("  Atlas: %d genes x %d columns\n", nrow(atlas), ncol(atlas)))

# Merge community assignments with atlas
merged <- merge(comm, atlas, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
cat(sprintf("  Merged: %d genes\n", nrow(merged)))

# Background gene universe = all genes in the network
universe <- comm$gene

# ======================================================================
# 2. Load MSigDB Hallmark gene sets
# ======================================================================
cat("\n--- Loading MSigDB Hallmark pathways ---\n")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("  Hallmark: %d gene sets\n", length(hallmark_sets)))

# ======================================================================
# 3. ORA helper function (Fisher's exact test)
# ======================================================================
run_ora <- function(target_genes, gene_set, universe) {
  # target_genes: genes in this community
  # gene_set: genes in the pathway/category
  # universe: all network genes
  target <- intersect(target_genes, universe)
  gs <- intersect(gene_set, universe)
  bg <- universe

  a <- length(intersect(target, gs))       # in community AND in set

b <- length(setdiff(target, gs))          # in community, NOT in set
  c <- length(setdiff(gs, target))          # NOT in community, in set
  d <- length(setdiff(bg, union(target, gs)))  # neither

  if (a == 0) {
    return(data.table(overlap = 0L, odds_ratio = 0, pvalue = 1))
  }

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "greater")

  data.table(
    overlap = as.integer(a),
    odds_ratio = as.numeric(ft$estimate),
    pvalue = as.numeric(ft$p.value)
  )
}

# ======================================================================
# 4. Run enrichment for a single community
# ======================================================================
enrich_community <- function(community_genes, community_id, level,
                             hallmark_sets, merged_dt, universe) {
  results <- list()
  n_genes <- length(community_genes)

  # -- A. Hallmark pathway enrichment --
  for (gs_name in names(hallmark_sets)) {
    res <- run_ora(community_genes, hallmark_sets[[gs_name]], universe)
    if (res$overlap > 0) {
      results[[length(results) + 1]] <- data.table(
        community_id = community_id,
        level = level,
        category = "hallmark",
        term = gs_name,
        overlap = res$overlap,
        n_genes = n_genes,
        odds_ratio = res$odds_ratio,
        pvalue = res$pvalue
      )
    }
  }

  # -- B. Drug target enrichment --
  drug_genes <- merged_dt[dgidb_druggable == TRUE | dgidb_druggable == "TRUE",
                          gene]
  drug_genes <- intersect(drug_genes, universe)
  if (length(drug_genes) > 0) {
    res <- run_ora(community_genes, drug_genes, universe)
    results[[length(results) + 1]] <- data.table(
      community_id = community_id, level = level,
      category = "drug_target", term = "DGIdb_druggable",
      overlap = res$overlap, n_genes = n_genes,
      odds_ratio = res$odds_ratio, pvalue = res$pvalue
    )
  }

  # -- C. COLOC gene enrichment (PP.H4 > 0.5) --
  if ("coloc_susie_best_pp4" %in% names(merged_dt)) {
    coloc_genes <- merged_dt[!is.na(coloc_susie_best_pp4) &
                               coloc_susie_best_pp4 > 0.5, gene]
    coloc_genes <- intersect(coloc_genes, universe)
    if (length(coloc_genes) > 0) {
      res <- run_ora(community_genes, coloc_genes, universe)
      results[[length(results) + 1]] <- data.table(
        community_id = community_id, level = level,
        category = "coloc", term = "COLOC_PP4_gt_0.5",
        overlap = res$overlap, n_genes = n_genes,
        odds_ratio = res$odds_ratio, pvalue = res$pvalue
      )
    }
  }

  # -- D. Sex class distribution (chi-squared vs background) --
  if ("sex_class" %in% names(merged_dt)) {
    comm_sex <- merged_dt[gene %in% community_genes & !is.na(sex_class) &
                            sex_class != "", .(gene, sex_class)]
    bg_sex   <- merged_dt[gene %in% universe & !is.na(sex_class) &
                            sex_class != "", .(gene, sex_class)]

    for (sx in c("Female_biased", "Male_biased", "Divergent")) {
      sx_genes <- bg_sex[sex_class == sx, gene]
      if (length(sx_genes) > 5) {
        res <- run_ora(community_genes, sx_genes, universe)
        results[[length(results) + 1]] <- data.table(
          community_id = community_id, level = level,
          category = "sex_class", term = sx,
          overlap = res$overlap, n_genes = n_genes,
          odds_ratio = res$odds_ratio, pvalue = res$pvalue
        )
      }
    }
  }

  # -- E. Progression genes --
  if ("has_progression_coloc" %in% names(merged_dt)) {
    prog_genes <- merged_dt[has_progression_coloc == TRUE |
                              has_progression_coloc == "TRUE", gene]
    prog_genes <- intersect(prog_genes, universe)
    if (length(prog_genes) > 0) {
      res <- run_ora(community_genes, prog_genes, universe)
      results[[length(results) + 1]] <- data.table(
        community_id = community_id, level = level,
        category = "progression", term = "progression_coloc",
        overlap = res$overlap, n_genes = n_genes,
        odds_ratio = res$odds_ratio, pvalue = res$pvalue
      )
    }
  }

  # -- F. Zonation class --
  if ("zonation_class" %in% names(merged_dt)) {
    for (zone in c("pericentral", "periportal")) {
      zone_genes <- merged_dt[zonation_class == zone, gene]
      zone_genes <- intersect(zone_genes, universe)
      if (length(zone_genes) > 5) {
        res <- run_ora(community_genes, zone_genes, universe)
        results[[length(results) + 1]] <- data.table(
          community_id = community_id, level = level,
          category = "zonation", term = zone,
          overlap = res$overlap, n_genes = n_genes,
          odds_ratio = res$odds_ratio, pvalue = res$pvalue
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  rbindlist(results)
}

# ======================================================================
# 5. Run enrichment across all macro and meso communities
# ======================================================================
cat("\n--- Running enrichment ---\n")

all_results <- list()
result_idx <- 0

for (level in c("macro", "meso")) {
  id_col <- paste0(level, "_id")
  community_ids <- sort(unique(merged[[id_col]]))
  cat(sprintf("  Level: %s (%d communities)\n", level, length(community_ids)))

  for (cid in community_ids) {
    comm_genes <- merged[get(id_col) == cid, gene]

    # Skip tiny communities (< 5 genes)
    if (length(comm_genes) < 5) next

    res <- enrich_community(
      community_genes = comm_genes,
      community_id = cid,
      level = level,
      hallmark_sets = hallmark_sets,
      merged_dt = merged,
      universe = universe
    )
    if (!is.null(res)) {
      result_idx <- result_idx + 1
      all_results[[result_idx]] <- res
    }
  }
}

if (length(all_results) == 0) {
  cat("WARNING: No enrichment results found. Check inputs.\n")
  enrichment <- data.table(
    community_id = integer(), level = character(), category = character(),
    term = character(), overlap = integer(), n_genes = integer(),
    odds_ratio = numeric(), pvalue = numeric(), padj = numeric()
  )
} else {
  enrichment <- rbindlist(all_results)

  # Multiple testing correction (BH across all tests)
  enrichment[, padj := p.adjust(pvalue, method = "BH")]

  # Sort by significance
  setorder(enrichment, padj, -odds_ratio)

  cat(sprintf("  Total enrichment tests: %d\n", nrow(enrichment)))
  cat(sprintf("  Significant (padj < 0.05): %d\n", sum(enrichment$padj < 0.05)))
}

# Save enrichment results
enrich_out <- file.path(COMM_DIR, "community_enrichment.csv")
fwrite(enrichment, enrich_out)
cat(sprintf("  Saved: %s\n", enrich_out))

# ======================================================================
# 6. Assign biological labels based on top enriched pathway
# ======================================================================
cat("\n--- Assigning community labels ---\n")

# Pathway-to-label mapping for cleaner names
pathway_labels <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"  = "Inflammatory-NF-kB",
  "HALLMARK_INFLAMMATORY_RESPONSE"      = "Inflammatory-Response",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"  = "Inflammatory-IFNg",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"  = "Inflammatory-IFNa",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"    = "Inflammatory-IL6",
  "HALLMARK_COMPLEMENT"                 = "Inflammatory-Complement",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"  = "Metabolic-OxPhos",
  "HALLMARK_FATTY_ACID_METABOLISM"      = "Metabolic-FattyAcid",
  "HALLMARK_BILE_ACID_METABOLISM"       = "Metabolic-BileAcid",
  "HALLMARK_XENOBIOTIC_METABOLISM"      = "Metabolic-Xenobiotic",
  "HALLMARK_ADIPOGENESIS"               = "Metabolic-Adipogenesis",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS"    = "Metabolic-Cholesterol",
  "HALLMARK_GLYCOLYSIS"                 = "Metabolic-Glycolysis",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Fibrotic-EMT",
  "HALLMARK_TGF_BETA_SIGNALING"        = "Fibrotic-TGFb",
  "HALLMARK_COAGULATION"                = "Fibrotic-Coagulation",
  "HALLMARK_APOPTOSIS"                  = "Stress-Apoptosis",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"  = "Stress-UPR",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "Stress-ROS",
  "HALLMARK_HYPOXIA"                    = "Stress-Hypoxia",
  "HALLMARK_P53_PATHWAY"                = "Stress-p53",
  "HALLMARK_DNA_REPAIR"                 = "Stress-DNARepair",
  "HALLMARK_E2F_TARGETS"               = "Proliferative-E2F",
  "HALLMARK_G2M_CHECKPOINT"             = "Proliferative-G2M",
  "HALLMARK_MYC_TARGETS_V1"            = "Proliferative-MYC",
  "HALLMARK_ANGIOGENESIS"               = "Vascular-Angiogenesis",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING" = "Signaling-WNT",
  "HALLMARK_NOTCH_SIGNALING"            = "Signaling-Notch",
  "HALLMARK_HEDGEHOG_SIGNALING"         = "Signaling-Hedgehog",
  "HALLMARK_KRAS_SIGNALING_UP"          = "Signaling-KRAS"
)

label_rows <- list()

for (level in c("macro", "meso")) {
  id_col <- paste0(level, "_id")
  community_ids <- sort(unique(merged[[id_col]]))

  for (cid in community_ids) {
    n_genes <- sum(merged[[id_col]] == cid)

    # Get hallmark enrichments for this community
    sub <- enrichment[level == ..level & community_id == cid &
                        category == "hallmark" & padj < 0.05]

    if (nrow(sub) > 0) {
      # Pick top pathway by odds ratio (among significant)
      setorder(sub, -odds_ratio)
      top_pathway <- sub$term[1]
      top_or <- sub$odds_ratio[1]
      top_padj <- sub$padj[1]

      # Map to clean label
      label <- pathway_labels[top_pathway]
      if (is.na(label)) {
        # Fallback: clean the Hallmark name
        label <- gsub("^HALLMARK_", "", top_pathway)
        label <- gsub("_", "-", label)
        label <- paste0(toupper(substring(label, 1, 1)),
                        tolower(substring(label, 2)))
      }
    } else {
      top_pathway <- "none"
      top_or <- NA_real_
      top_padj <- NA_real_
      label <- "Uncharacterized"
    }

    label_rows[[length(label_rows) + 1]] <- data.table(
      community_id = cid,
      level = level,
      label = label,
      top_pathway = top_pathway,
      top_odds_ratio = top_or,
      top_padj = top_padj,
      n_genes = n_genes
    )
  }
}

labels_dt <- rbindlist(label_rows)
labels_out <- file.path(COMM_DIR, "community_labels.csv")
fwrite(labels_dt, labels_out)
cat(sprintf("  Labels saved: %s (%d rows)\n", labels_out, nrow(labels_dt)))

# Print label summary
cat("\n  Macro labels:\n")
macro_labels <- labels_dt[level == "macro"]
setorder(macro_labels, community_id)
for (i in seq_len(nrow(macro_labels))) {
  r <- macro_labels[i]
  cat(sprintf("    M%d: %s (%d genes, top: %s, OR=%.1f)\n",
              r$community_id, r$label, r$n_genes, r$top_pathway,
              ifelse(is.na(r$top_odds_ratio), 0, r$top_odds_ratio)))
}

# ======================================================================
# 7. Build structured JSON for portal
# ======================================================================
cat("\n--- Building communities JSON ---\n")

build_level_json <- function(level, merged_dt, labels_dt, enrichment_dt) {
  id_col <- paste0(level, "_id")
  community_ids <- sort(unique(merged_dt[[id_col]]))
  level_list <- list()

  for (cid in community_ids) {
    genes <- merged_dt[get(id_col) == cid, gene]

    lab_row <- labels_dt[level == ..level & community_id == cid]
    label <- if (nrow(lab_row) > 0) lab_row$label[1] else "Uncharacterized"
    top_pw <- if (nrow(lab_row) > 0) lab_row$top_pathway[1] else "none"

    # Top 10 enrichment results for this community
    enrich_sub <- enrichment_dt[level == ..level & community_id == cid & padj < 0.05]
    setorder(enrich_sub, padj)
    enrich_top <- head(enrich_sub, 10)

    enrich_list <- lapply(seq_len(nrow(enrich_top)), function(j) {
      list(
        category = enrich_top$category[j],
        term = enrich_top$term[j],
        odds_ratio = round(enrich_top$odds_ratio[j], 2),
        padj = signif(enrich_top$padj[j], 3)
      )
    })

    level_list[[length(level_list) + 1]] <- list(
      id = cid,
      label = label,
      n_genes = length(genes),
      genes = as.list(genes),
      top_pathway = top_pw,
      enrichment = enrich_list
    )
  }

  level_list
}

portal_json <- list(
  macro = build_level_json("macro", merged, labels_dt, enrichment),
  meso  = build_level_json("meso",  merged, labels_dt, enrichment),
  micro = NULL  # micro communities listed in assignments CSV only (too many for JSON)
)

# Micro summary: just id, n_genes, parent meso/macro
micro_ids <- sort(unique(merged$micro_id))
micro_summary <- lapply(micro_ids, function(mid) {
  genes <- merged[micro_id == mid, gene]
  meso_parent <- merged[micro_id == mid, meso_id][1]
  macro_parent <- merged[micro_id == mid, macro_id][1]
  list(
    id = mid,
    n_genes = length(genes),
    meso_parent = meso_parent,
    macro_parent = macro_parent
  )
})
portal_json$micro <- micro_summary

json_out <- file.path(COMM_DIR, "communities.json")
write_json(portal_json, json_out, auto_unbox = TRUE, pretty = FALSE)
cat(sprintf("  JSON saved: %s\n", json_out))

# ======================================================================
# Summary
# ======================================================================
cat(sprintf("\n=== 268 Complete ===\n"))
cat(sprintf("  Enrichment results: %d tests, %d significant (padj < 0.05)\n",
            nrow(enrichment),
            sum(enrichment$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Labels: %d macro + %d meso communities labeled\n",
            sum(labels_dt$level == "macro"),
            sum(labels_dt$level == "meso")))
cat(sprintf("  Outputs:\n"))
cat(sprintf("    %s\n", enrich_out))
cat(sprintf("    %s\n", labels_out))
cat(sprintf("    %s\n", json_out))
cat(sprintf("  End: %s\n", format(Sys.time())))
