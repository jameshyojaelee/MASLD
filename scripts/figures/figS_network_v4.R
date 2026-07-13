#!/usr/bin/env Rscript
# ============================================================================
# figS_network_v4.R — consolidated supplementary panels for v4 Architecture C
# ----------------------------------------------------------------------------
# Produces 7 PDFs in figures/supplementary/figS_network/:
#   A architecture schematic
#   B per-community edge-type prevalence heatmap (top-20 F3-F4)
#   C LOCO replication distribution across D-F2 edges
#   D honest STRING benchmark (Jaccard + per-layer fold enrichment)
#   E mid-stage (F1-F3) inflection hub subnetworks (emerging vs dissolving)
#   F gene-neighborhood exemplars (THRB, HNF4A, TM6SF2, CFLAR, HSD17B13,
#     SERPINA1)
#   G portal query summaries (THRB, Q2 COLOC ∩ mid-stage (F1-F3) inflection, Q3 fibrogenic
#     druggable)
#
# Data sources: see CLAUDE.md and memory/network_final_decisions.md.
# ============================================================================

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(purrr)
    library(ggplot2); library(ggrepel); library(patchwork); library(jsonlite)
    library(igraph); library(ggraph); library(tidygraph)
    library(scales); library(stringr); library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

NET_DIR     <- file.path(BASE, "RNA-seq/results/network")
COMM_DIR    <- file.path(NET_DIR, "communities_f2")
PORTAL_DIR  <- file.path(NET_DIR, "portal_export_v2")
GRAPH_DIR   <- file.path(PORTAL_DIR, "gene_graphs")
BUNDLE_DIR  <- file.path(PORTAL_DIR, "query_bundles")
BENCH_DIR   <- file.path(BASE, "figures/misc/network_benchmark_v2/data")
OUT_DIR     <- FIGS_NET_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

edge_type_colors <- c(
    S            = "#9E9E9E",
    `D-F2`       = "#EF6C00",
    `D-COLOC`    = "#6A1B9A",
    `D-LR`       = "#2E7D32",
    `D-ceRNA`    = "#00838F",
    `D-XSpecies` = "#AD1457"
)

# ============================================================================
# Panel A — architecture schematic
# ============================================================================
panel_A <- function() {
    nodes <- tibble(
        id   = c("STRING", "D-F2", "D-COLOC", "D-LR", "D-ceRNA", "D-XSpecies"),
        kind = c("backbone", rep("dedge", 5)),
        label = c(
            # Edge counts from the current network build (290_edge_annotation_atlas.py,
            # edge_annotation_summary.csv: 118,791 STRING >=700 edges total; "contested"
            # = STRING edge that also carries a disease-native D edge). The previous
            # "115,020 / 3,626" pair was stale (3,626 coincided with an old DEG count).
            "STRING v12 >=700\nbackbone\n104,389 S-only\n14,402 contested",
            "D-F2\nstage-stratified\ndelta-r coexpression\n3,260,461 edges",
            "D-COLOC\nGWAS locus pairs\n111 edges\n(annotation)",
            "D-LR\nLIANA F2-differential\n9 edges",
            "D-ceRNA\nMASLD lncRNA\n0 edges (DB sparsity)",
            "D-XSpecies\nFPC switch\n(mouse-human)"
        ),
        x = c(0,  2.4, 1.8, -1.8, -2.4,  0.0),
        y = c(0,  0.9, -1.5, -1.5,  0.9,  2.2)
    )
    segs <- tibble(
        xend = 0, yend = 0,
        x    = nodes$x[-1], y = nodes$y[-1],
        type = nodes$id[-1]
    )
    ann <- tibble(
        x = -3.3, y = -3.2,
        label = paste(
            "Annotation atlas (3,375,599 edges x 40+ columns):",
            "F-stage, cell-type driver, cohort support, COLOC PP4,",
            "druggability, sex, conserved_mouse, cascade_module",
            sep = "\n")
    )
    p <- ggplot() +
        geom_segment(data = segs,
                     aes(x = x, y = y, xend = xend, yend = yend, color = type),
                     linewidth = 1.4, lineend = "round") +
        geom_point(data = nodes,
                   aes(x = x, y = y, fill = kind),
                   shape = 21, size = 16, color = "grey20", stroke = 0.8) +
        geom_text(data = nodes, aes(x = x, y = y, label = label),
                  size = GEOM_TEXT_6PT, lineheight = 0.95) +
        geom_text(data = ann, aes(x = x, y = y, label = label),
                  hjust = 0, size = GEOM_TEXT_6PT, color = "black") +
        scale_color_manual(values = edge_type_colors, name = "Edge type") +
        scale_fill_manual(values = c(backbone = "#E0E0E0", dedge = "#FFF3E0"),
                          guide = "none") +
        coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3)) +
        theme_masld() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank())
    message("[caption] Architecture C: STRING backbone + 5 disease-native edge types. ",
            "Type-D edges carry biology STRING cannot have.")
    ggsave(file.path(OUT_DIR, "figS_network_A_architecture.pdf"), p,
           width = 180/25.4, height = 120/25.4, device = cairo_pdf)
    message("Wrote figS_network_A_architecture.pdf")
}

# ============================================================================
# Panel B — edge-type prevalence heatmap across top-20 F3-F4 communities
# ============================================================================
panel_B <- function() {
    message("Panel B: loading edge atlas (parquet)...")
    atlas <- data.table::fread(file.path(NET_DIR, "edge_annotation_atlas.csv"))
    f34 <- read_csv(file.path(COMM_DIR, "communities_F34.csv"),
                    show_col_types = FALSE)
    labs_json <- fromJSON(file.path(COMM_DIR, "community_labels.json"),
                          simplifyVector = FALSE)

    type_col <- NA_character_
    for (cand in c("edge_type", "type", "layer")) {
        if (cand %in% names(atlas)) { type_col <- cand; break }
    }
    if (is.na(type_col)) {
        # derive from flags
        flag_map <- list(
            `D-F2`       = grep("^is_d_f2$|^d_f2$",       names(atlas), value = TRUE)[1],
            `D-COLOC`    = grep("^is_d_coloc$|^d_coloc$", names(atlas), value = TRUE)[1],
            `D-LR`       = grep("^is_d_lr$|^d_lr$",       names(atlas), value = TRUE)[1],
            `D-ceRNA`    = grep("^is_d_cerna$|^d_cerna$", names(atlas), value = TRUE)[1],
            `D-XSpecies` = grep("^is_d_xspecies$|^d_xs$", names(atlas), value = TRUE)[1],
            S            = grep("^is_string$|^in_string", names(atlas), value = TRUE)[1]
        )
        flag_map <- flag_map[!is.na(flag_map)]
        if (length(flag_map) == 0) {
            stop("Panel B: cannot find edge-type columns in atlas. ",
                 "Columns available: ",
                 paste(head(names(atlas), 20), collapse = ", "))
        }
        message("  deriving edge_type from flag columns: ",
                paste(names(flag_map), collapse = ", "))
        # priority D-F2 > D-COLOC > D-LR > D-ceRNA > D-XSpecies > S
        order_types <- intersect(
            c("D-F2","D-COLOC","D-LR","D-ceRNA","D-XSpecies","S"),
            names(flag_map))
        edge_type <- rep("S", nrow(atlas))
        for (t in rev(order_types)) {
            col <- flag_map[[t]]
            v <- atlas[[col]]
            idx <- which(!is.na(v) & as.logical(v))
            edge_type[idx] <- t
        }
        atlas$edge_type <- edge_type
        type_col <- "edge_type"
    }

    top20 <- f34 |> count(community_id, name = "size") |>
        arrange(desc(size)) |> head(20)
    f34_map <- f34 |> filter(community_id %in% top20$community_id) |>
        select(gene, community_id)

    a_small <- atlas |>
        select(gene_a, gene_b, all_of(type_col)) |>
        inner_join(f34_map |> rename(gene_a = gene, cid_a = community_id),
                   by = "gene_a") |>
        inner_join(f34_map |> rename(gene_b = gene, cid_b = community_id),
                   by = "gene_b") |>
        filter(cid_a == cid_b) |>
        transmute(community_id = cid_a,
                  edge_type = .data[[type_col]])

    prev <- a_small |>
        count(community_id, edge_type, name = "n_edges") |>
        group_by(community_id) |>
        mutate(frac = n_edges / sum(n_edges)) |>
        ungroup()

    lab_df <- purrr::map_dfr(names(labs_json$F34), function(cid) {
        tibble(community_id = as.integer(cid),
               label = labs_json$F34[[cid]]$label %||% "Mixed",
               n_genes = labs_json$F34[[cid]]$n_genes %||% NA_integer_)
    })
    prev <- prev |> left_join(lab_df, by = "community_id") |>
        mutate(row = sprintf("C%d %s (n=%d)", community_id, label, n_genes))

    ord <- top20 |> left_join(lab_df, by = "community_id") |>
        mutate(row = sprintf("C%d %s (n=%d)", community_id, label, n_genes)) |>
        arrange(desc(size)) |> pull(row)
    prev$row <- factor(prev$row, levels = rev(ord))

    p <- ggplot(prev, aes(x = edge_type, y = row, fill = frac)) +
        geom_tile(color = "white", linewidth = 0.3) +
        geom_text(aes(label = scales::percent(frac, accuracy = 0.1)),
                  size = GEOM_TEXT_6PT, color = "black") +
        scale_fill_gradient(low = "#F5F5F5", high = "#C2185B",
                            labels = percent_format(accuracy = 1),
                            name = "Fraction of\ninternal edges") +
        labs(x = "Edge type", y = NULL) +
        theme_masld() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    message("[caption] Per-community edge-type prevalence (F3-F4 top 20)")
    ggsave(file.path(OUT_DIR, "figS_network_B_evidence_diversity.pdf"), p,
           width = 180/25.4, height = 160/25.4, device = cairo_pdf)
    message("Wrote figS_network_B_evidence_diversity.pdf")
}

# ============================================================================
# Panel C — LOCO replication
# ============================================================================
panel_C <- function() {
    loco <- read_csv(file.path(NET_DIR, "edges_d_f2_loco.csv"),
                     show_col_types = FALSE)
    frac <- loco$replicates_loco_fraction
    frac <- frac[is.finite(frac)]
    med  <- median(frac)
    pass <- mean(frac >= 0.70)

    p <- ggplot(tibble(f = frac), aes(x = f)) +
        geom_histogram(bins = 60, fill = "#FFCC80", color = "grey25",
                       linewidth = 0.2) +
        geom_vline(xintercept = med, color = "#C2185B",
                   linewidth = 0.8, linetype = "dashed") +
        geom_vline(xintercept = 0.70, color = "#1565C0",
                   linewidth = 0.8, linetype = "dotted") +
        annotate("text", x = med, y = Inf,
                 label = sprintf("median = %.3f", med),
                 vjust = 2, hjust = -0.05, color = "#C2185B", size = GEOM_TEXT_6PT) +
        annotate("text", x = 0.70, y = Inf,
                 label = sprintf("threshold = 0.70\n%.1f%% pass",
                                 100 * pass),
                 vjust = 2, hjust = 1.05, color = "#1565C0", size = GEOM_TEXT_6PT) +
        labs(x = "Fraction of LOCO folds replicated",
             y = "Number of edges") +
        theme_masld()
    message(sprintf("[caption] LOCO replication of D-F2 edges (%s edges across leave-one-cohort-out folds)",
                    format(length(frac), big.mark = ",")))
    ggsave(file.path(OUT_DIR, "figS_network_C_loco_replication.pdf"), p,
           width = 140/25.4, height = 95/25.4, device = cairo_pdf)
    message("Wrote figS_network_C_loco_replication.pdf")
}

# ============================================================================
# Panel D — honest benchmark
# ============================================================================
panel_D <- function() {
    jac_df <- tibble(
        comparison = "D-F2 vs STRING >= 700",
        jaccard    = 0.0011
    )
    p_jac <- ggplot(jac_df, aes(x = comparison, y = jaccard)) +
        geom_col(fill = "#EF6C00", width = 0.45) +
        geom_hline(yintercept = 0.30, linetype = "dashed",
                   color = "#1565C0", linewidth = 0.7) +
        annotate("text", x = 1, y = 0.32,
                 label = "Threshold = 0.30 (pre-reg)",
                 color = "#1565C0", size = GEOM_TEXT_6PT, hjust = 0.5) +
        annotate("text", x = 1, y = 0.02,
                 label = "0.0011 — D-F2 captures\nbiology STRING cannot",
                 size = GEOM_TEXT_6PT, hjust = 0.5, color = "black") +
        scale_y_continuous(limits = c(0, 0.35), expand = expansion(0)) +
        labs(x = NULL, y = "Jaccard") +
        theme_masld() +
        theme(axis.text.x = element_text(size = 6))

    bench <- read_csv(file.path(BENCH_DIR, "bench_per_layer.csv"),
                      show_col_types = FALSE)
    keep_sets <- c("clinical_drug", "govaere25",
                   "resmetirom_pathway", "resmetirom")
    bench2 <- bench |> filter(gold_set %in% keep_sets)
    bench2$gold_set <- factor(bench2$gold_set,
        levels = intersect(keep_sets, unique(bench2$gold_set)))

    p_bench <- ggplot(bench2,
        aes(x = layer, y = fold_over_null, fill = gold_set)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.75,
                 color = "grey25", linewidth = 0.2) +
        geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
        scale_fill_brewer(palette = "Set2", name = "Gold set") +
        labs(x = "Edge layer",
             y = "Fold over degree-preserving null") +
        theme_masld() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              legend.position = "bottom")

    combined <- (p_jac | p_bench) + plot_layout(widths = c(0.6, 1.4))
    message("[caption] (i) D-F2 vs STRING overlap. ",
            "(ii) Per-layer fold-enrichment on MASLD gold sets — ",
            "STRING wins generic; D-layers win disease-native.")
    ggsave(file.path(OUT_DIR, "figS_network_D_honest_benchmark.pdf"),
           combined, width = 180/25.4, height = 110/25.4, device = cairo_pdf)
    message("Wrote figS_network_D_honest_benchmark.pdf")
}

# ============================================================================
# Panel E — mid-stage (F1-F3) inflection hub subnetworks
# ============================================================================
panel_E <- function() {
    e <- read_csv(file.path(NET_DIR, "edges_d_f2_loco.csv"),
                  show_col_types = FALSE) |>
        select(gene_a, gene_b, emergence_stage, replicates_loco_fraction)

    build_hub_graph <- function(edges, stage, top_k = 20) {
        ed <- edges |>
            filter(emergence_stage == stage,
                   is.na(replicates_loco_fraction) |
                       replicates_loco_fraction >= 0.70)
        if (nrow(ed) == 0) return(NULL)
        deg <- bind_rows(
            ed |> count(gene = gene_a),
            ed |> count(gene = gene_b)
        ) |> group_by(gene) |> summarise(n = sum(n), .groups = "drop") |>
            arrange(desc(n)) |> head(top_k)
        keep <- deg$gene
        g_ed <- ed |> filter(gene_a %in% keep, gene_b %in% keep) |>
            select(from = gene_a, to = gene_b) |> distinct()
        if (nrow(g_ed) < 5) {
            g_ed <- ed |> filter(gene_a %in% keep | gene_b %in% keep) |>
                select(from = gene_a, to = gene_b) |> distinct() |>
                head(200)
        }
        nodes <- tibble(name = unique(c(g_ed$from, g_ed$to))) |>
            left_join(deg |> rename(name = gene), by = "name") |>
            mutate(n = dplyr::coalesce(n, 1),
                   hub = name %in% keep)
        tidygraph::tbl_graph(nodes = nodes, edges = g_ed, directed = FALSE)
    }

    g_em  <- build_hub_graph(e, "F2_emerging")
    g_dis <- build_hub_graph(e, "F2_dissolving")

    render <- function(g, title_txt) {
        if (is.null(g)) {
            return(ggplot() + theme_void())
        }
        ggraph(g, layout = "fr") +
            geom_edge_link(alpha = 0.3, edge_width = 0.3, color = "#EF6C00") +
            geom_node_point(aes(size = n, filter = hub),
                            fill = "#C2185B", shape = 21, color = "grey20") +
            geom_node_point(aes(size = n, filter = !hub),
                            fill = "#E0E0E0", shape = 21, color = "grey30",
                            alpha = 0.7) +
            geom_node_text(aes(label = name, filter = hub),
                           size = GEOM_TEXT_6PT, repel = TRUE, color = "black") +
            scale_size(range = c(2, 7), guide = "none") +
            theme_masld() +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank(), panel.grid = element_blank())
    }

    p1 <- render(g_em,  "F2-emerging hubs (top 20)")
    p2 <- render(g_dis, "F2-dissolving hubs (top 20)")
    combined <- p1 | p2
    message("[caption] Left: F2-emerging hubs (top 20). Right: F2-dissolving hubs (top 20).")
    ggsave(file.path(OUT_DIR, "figS_network_E_switch_hubs.pdf"),
           combined, width = 180/25.4, height = 110/25.4, device = cairo_pdf)
    message("Wrote figS_network_E_switch_hubs.pdf")
}

# ============================================================================
# Panel F — 6 gene neighborhoods
# ============================================================================
panel_F <- function() {
    exemplars <- c("THRB", "HNF4A", "TM6SF2", "CFLAR", "HSD17B13", "SERPINA1")
    build_one <- function(sym) {
        fp <- file.path(GRAPH_DIR, paste0(sym, ".json"))
        if (!file.exists(fp)) {
            message("  missing gene graph for ", sym); return(NULL)
        }
        d <- fromJSON(fp, simplifyVector = FALSE)
        collect <- list()
        for (subset_name in names(d$neighbors)) {
            lst <- d$neighbors[[subset_name]]
            if (length(lst) == 0) next
            type_key <- dplyr::case_when(
                str_detect(subset_name, "string")  ~ "S",
                str_detect(subset_name, "d_f2")    ~ "D-F2",
                str_detect(subset_name, "d_coloc") ~ "D-COLOC",
                str_detect(subset_name, "d_lr")    ~ "D-LR",
                str_detect(subset_name, "d_cerna") ~ "D-ceRNA",
                str_detect(subset_name, "d_x")     ~ "D-XSpecies",
                TRUE ~ "S"
            )
            for (n in lst) {
                collect[[length(collect) + 1]] <- tibble(
                    partner = n$partner %||% NA_character_,
                    type    = type_key,
                    score   = n$string_score %||% NA_real_
                )
            }
        }
        if (length(collect) == 0) return(NULL)
        nb <- bind_rows(collect) |> filter(!is.na(partner)) |>
            group_by(partner) |>
            slice_max(order_by = dplyr::coalesce(score, 0.9),
                      n = 1, with_ties = FALSE) |>
            ungroup() |>
            mutate(rank_score = dplyr::coalesce(score, 0.9)) |>
            arrange(desc(rank_score)) |> head(25)
        edges <- tibble(from = sym, to = nb$partner, type = nb$type)
        nodes <- tibble(name = unique(c(sym, nb$partner))) |>
            mutate(hub = name == sym)
        tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
    }

    panels <- purrr::map(exemplars, function(sym) {
        g <- build_one(sym)
        if (is.null(g)) return(ggplot() + theme_void())
        ggraph(g, layout = "fr") +
            geom_edge_link(aes(color = type), edge_width = 0.5, alpha = 0.8) +
            geom_node_point(aes(fill = hub), shape = 21, size = 3.5,
                            color = "grey25") +
            geom_node_text(aes(label = name), repel = TRUE, size = GEOM_TEXT_6PT,
                           color = "black") +
            scale_edge_color_manual(values = edge_type_colors,
                                    name = "Edge type", drop = FALSE) +
            scale_fill_manual(values = c(`TRUE` = "#C2185B",
                                         `FALSE` = "#ECEFF1"),
                              guide = "none") +
            theme_masld() +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank(), panel.grid = element_blank(),
                  legend.position = "bottom")
    })
    combined <- wrap_plots(panels, ncol = 3) +
        plot_layout(guides = "collect") +
        plot_annotation(theme = theme(legend.position = "bottom"))
    message(sprintf("[caption] Gene-neighborhood exemplars (top-25 partners): %s",
                    paste(exemplars, collapse = ", ")))
    ggsave(file.path(OUT_DIR, "figS_network_F_neighborhood_exemplars.pdf"),
           combined, width = 180/25.4, height = 180/25.4, device = cairo_pdf)
    message("Wrote figS_network_F_neighborhood_exemplars.pdf")
}

# ============================================================================
# Panel G — portal queries
# ============================================================================
panel_G <- function() {
    thrb_graph <- fromJSON(file.path(GRAPH_DIR, "THRB.json"),
                           simplifyVector = FALSE)
    q2 <- fromJSON(file.path(BUNDLE_DIR, "Q2_coloc_f2_switch.json"),
                   simplifyVector = FALSE)
    q3 <- fromJSON(file.path(BUNDLE_DIR, "Q3_fibrogenic_druggable.json"),
                   simplifyVector = FALSE)

    # Q1 — THRB STRING neighbors (top 10 by score)
    ns <- thrb_graph$neighbors$neighbors_string %||% list()
    nb_tbl <- map_dfr(ns, function(x) tibble(
        partner = x$partner %||% NA_character_,
        score   = x$string_score %||% NA_real_))
    ec <- thrb_graph$edge_counts %||% list()
    n_total <- sum(unlist(ec))
    top10 <- nb_tbl |> filter(!is.na(partner)) |>
        arrange(desc(score)) |> head(10)
    p_thrb <- ggplot(top10,
        aes(x = reorder(partner, score), y = score)) +
        geom_col(fill = "#C2185B") + coord_flip() +
        labs(x = NULL, y = "STRING score") +
        theme_masld()

    # Q2/Q3 genes are character vectors with `n` count
    mk_panel <- function(q, color, label) {
        genes <- unlist(q$genes) %||% character(0)
        n_tot <- q$n %||% length(genes)
        top <- head(genes, 10)
        df <- tibble(gene = top, rank = seq_along(top))
        ggplot(df, aes(x = reorder(gene, -rank), y = 1)) +
            geom_col(fill = color) + coord_flip() +
            labs(x = NULL, y = NULL) +
            theme_masld() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
    }
    p_q2 <- mk_panel(q2, "#6A1B9A", "Q2. COLOC \u2229 mid-stage (F1-F3) inflection")
    p_q3 <- mk_panel(q3, "#2E7D32", "Q3. Fibrogenic druggable")

    n2 <- q2$n %||% length(unlist(q2$genes) %||% character(0))
    n3 <- q3$n %||% length(unlist(q3$genes) %||% character(0))
    combined <- (p_thrb | p_q2 | p_q3)
    message(sprintf(
        "[caption] Portal-ready query bundles. Q1. THRB neighborhood (n=%d edges, top 10 STRING partners by score). ",
        n_total),
        sprintf("Q2. COLOC \u2229 mid-stage (F1-F3) inflection (n=%d genes, first 10 members). ", n2),
        sprintf("Q3. Fibrogenic druggable (n=%d genes, first 10 members).", n3))
    ggsave(file.path(OUT_DIR, "figS_network_G_portal_queries.pdf"),
           combined, width = 180/25.4, height = 95/25.4, device = cairo_pdf)
    message("Wrote figS_network_G_portal_queries.pdf")
}

# ============================================================================
# Main
# ============================================================================
main <- function() {
    wp <- Sys.getenv("FIGS_NET_PANELS", "ACDGFEB")
    run_panel <- function(letter, fn) {
        if (!grepl(letter, wp)) return(invisible())
        tryCatch(fn(), error = function(e)
            message("Panel ", letter, " FAILED: ", conditionMessage(e)))
    }
    run_panel("A", panel_A)
    run_panel("C", panel_C)
    run_panel("D", panel_D)
    run_panel("G", panel_G)
    run_panel("F", panel_F)
    run_panel("E", panel_E)
    run_panel("B", panel_B)  # heaviest — parquet
    message("figS_network_v4 complete. Outputs in: ", OUT_DIR)
}

if (!interactive()) main()
