#!/usr/bin/env Rscript
# KEY MESSAGE: Gene-membership overlap organizes frozen pathways and Hotspot programs into outcome-blind molecular systems before stage or continuum results are read.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(irlba)
  library(igraph)
  library(uwot)
})

options(digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 66_freeze_molecular_systems_geometry.R ANALYSIS_CANDIDATE OUTPUT_CANDIDATE",
    call. = FALSE
  )
}
analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()
set.seed(contract$seed)

geometry_dir <- file.path(output_candidate, "geometry")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(geometry_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

family_order <- c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc")
family_labels <- c(
  hallmark = "Hallmark", kegg = "KEGG", reactome = "Reactome",
  go_bp = "GO biological process", go_mf = "GO molecular function",
  go_cc = "GO cellular component", hotspot = "Hotspot"
)
expected_sizes <- unlist(contract$pathway_collections[family_order], use.names = TRUE)
evaluation_cohorts <- contract$evaluation_cohorts

clean_set_label <- function(set_id, collection, description = NA_character_) {
  if (collection == "reactome" && is.character(description) &&
      length(description) == 1L && !is.na(description) && nzchar(description) &&
      nchar(description) <= 70L) {
    answer <- description
  } else {
    answer <- sub(
      "^(HALLMARK_|KEGG_|REACTOME_|GOBP_|GOMF_|GOCC_)", "", set_id
    )
    answer <- tools::toTitleCase(tolower(gsub("_", " ", answer, fixed = TRUE)))
  }
  answer <- gsub("Nfkb", "NF-kB", answer, fixed = TRUE)
  answer <- gsub("Tgf Beta", "TGF-beta", answer, fixed = TRUE)
  answer
}

collapse_testability <- function(path, ids, id_column, collection, expected_size) {
  x <- fread(path, na.strings = c("", "NA"))
  ml_assert(all(c(id_column, "dataset", "testable") %in% names(x)),
            paste0("Testability schema drift: ", collection))
  if (!"testability_reason" %in% names(x)) {
    x[, testability_reason := fifelse(
      testable %in% TRUE, "testable",
      "minimum_gene_or_retained_membership_weight_gate_not_met"
    )]
  }
  setnames(x, id_column, "feature_id")
  x <- x[dataset %in% evaluation_cohorts]
  ml_assert(nrow(x) == expected_size * length(evaluation_cohorts),
            paste0("Testability family size drift: ", collection))
  ml_assert(!anyDuplicated(x[, .(feature_id, dataset)]),
            paste0("Duplicated testability rows: ", collection))
  summary <- x[, .(
    n_evaluation_cohorts = uniqueN(dataset),
    n_testable_evaluation_cohorts = sum(testable %in% TRUE),
    geometry_eligible = uniqueN(dataset) == length(evaluation_cohorts) &&
      all(testable %in% TRUE),
    testability_reason = if (all(testable %in% TRUE)) "testable_in_both_evaluation_cohorts"
      else paste(sort(unique(testability_reason[testable %in% FALSE])), collapse = ";")
  ), by = feature_id]
  answer <- merge(
    data.table(feature_id = ids), summary, by = "feature_id", all.x = TRUE,
    sort = FALSE
  )
  answer[, `:=`(
    collection = collection,
    collection_label = family_labels[[collection]],
    complete_family_size = expected_size
  )]
  answer
}

pathway_features <- list()
pathway_memberships <- list()
testability_rows <- list()
input_rows <- list()

pathway_root <- ml_resolve(contract$pathway_dir)
for (collection in family_order) {
  collection_id <- collection
  family_size <- expected_sizes[[collection]]
  gene_path <- file.path(pathway_root, paste0(collection, "_term2gene.csv"))
  name_path <- file.path(pathway_root, paste0(collection, "_term2name.csv"))
  test_path <- file.path(
    analysis_candidate, "pathway_tf", "pathways", collection, "testability.tsv"
  )
  ml_assert(all(file.exists(c(gene_path, name_path, test_path))),
            paste0("Missing frozen pathway input: ", collection))
  genes <- fread(gene_path, na.strings = c("", "NA"))
  names_table <- fread(name_path, na.strings = c("", "NA"))
  setnames(genes, c("gs_name", "gene_symbol"), c("feature_id", "source_gene"),
           skip_absent = TRUE)
  setnames(names_table, c("gs_name", "gs_description"),
           c("feature_id", "description"), skip_absent = TRUE)
  ml_assert(all(c("feature_id", "source_gene") %in% names(genes)),
            paste0("TERM2GENE schema drift: ", collection))
  ids <- unique(genes$feature_id)
  ml_assert(length(ids) == family_size,
            paste0("TERM2GENE family size drift: ", collection))
  testability <- collapse_testability(
    test_path, ids, "set_id", collection, family_size
  )
  testability_rows[[collection]] <- testability
  eligible_ids <- testability[geometry_eligible == TRUE, feature_id]
  registry <- unique(genes[, .(feature_id)])
  registry <- merge(registry, names_table, by = "feature_id", all.x = TRUE, sort = FALSE)
  registry[, `:=`(
    collection = collection,
    collection_label = family_labels[[collection]],
    feature_type = "pathway",
    cell_type = NA_character_,
    module = NA_character_
  )]
  registry[, display_label := vapply(
    seq_len(.N),
    function(index) clean_set_label(feature_id[[index]], collection_id,
                                    description[[index]]),
    character(1)
  )]
  pathway_features[[collection]] <- registry[feature_id %in% eligible_ids]
  pathway_memberships[[collection]] <- unique(
    genes[feature_id %in% eligible_ids,
          .(feature_id, source_gene = toupper(trimws(source_gene)), collection)]
  )
  input_rows[[length(input_rows) + 1L]] <- data.table(
    path = normalizePath(c(gene_path, name_path, test_path), mustWork = TRUE),
    input_role = c("pathway_membership", "pathway_name", "outcome_blind_testability"),
    collection = collection
  )
}

hotspot_membership_path <- file.path(
  analysis_candidate, "programs", "hotspot",
  "hotspot_signature_excluded_membership.tsv.gz"
)
hotspot_testability_path <- file.path(
  analysis_candidate, "programs", "hotspot", "hotspot_testability.tsv"
)
ml_assert(all(file.exists(c(hotspot_membership_path, hotspot_testability_path))),
          "Missing frozen Hotspot geometry input")
hotspot_membership <- fread(hotspot_membership_path, na.strings = c("", "NA"))
required_hotspot <- c(
  "program_uid", "gene_symbol", "retained_relative_weight",
  "excluded_signature_gene", "cell_type", "module", "module_name"
)
ml_assert(all(required_hotspot %in% names(hotspot_membership)),
          "Hotspot membership schema drift")
hotspot_ids <- sort(unique(hotspot_membership$program_uid))
ml_assert(length(hotspot_ids) == contract$program_family_size,
          "Hotspot complete family size drift")
hotspot_testability <- collapse_testability(
  hotspot_testability_path, hotspot_ids, "program_uid", "hotspot",
  contract$program_family_size
)
testability_rows[["hotspot"]] <- hotspot_testability
eligible_hotspot <- hotspot_testability[geometry_eligible == TRUE, feature_id]
hotspot_features <- unique(hotspot_membership[
  program_uid %in% eligible_hotspot,
  .(
    feature_id = program_uid, collection = "hotspot",
    collection_label = family_labels[["hotspot"]], feature_type = "hotspot_program",
    display_label = module_name, description = module_name,
    cell_type, module = as.character(module)
  )
])
hotspot_memberships <- hotspot_membership[
  program_uid %in% eligible_hotspot & excluded_signature_gene == FALSE,
  .(
    feature_id = program_uid, source_gene = toupper(trimws(gene_symbol)),
    collection = "hotspot", source_weight = as.numeric(retained_relative_weight)
  )
]
input_rows[[length(input_rows) + 1L]] <- data.table(
  path = normalizePath(c(hotspot_membership_path, hotspot_testability_path),
                       mustWork = TRUE),
  input_role = c("hotspot_membership", "outcome_blind_testability"),
  collection = "hotspot"
)

feature_registry <- rbindlist(
  c(pathway_features, list(hotspot = hotspot_features)),
  use.names = TRUE, fill = TRUE
)
ml_assert(!anyDuplicated(feature_registry$feature_id),
          "Feature IDs collide across molecular-layer collections")
expected_eligible <- c(
  hallmark = 44L, kegg = 260L, reactome = 1027L, go_bp = 2811L,
  go_mf = 720L, go_cc = 467L, hotspot = 116L
)
eligible_counts <- feature_registry[, .N, by = collection]
ml_assert(
  identical(
    eligible_counts[match(names(expected_eligible), collection), N],
    unname(expected_eligible)
  ),
  "Outcome-blind eligible feature count drift"
)
ml_assert(nrow(feature_registry) == 5445L,
          "Expected 5,445 outcome-blind molecular-system nodes")

complete_ledger <- rbindlist(testability_rows, use.names = TRUE, fill = TRUE)
ml_assert(nrow(complete_ledger) ==
            sum(expected_sizes) + contract$program_family_size,
          "Complete molecular-system family ledger drift")
ml_assert(sum(complete_ledger$geometry_eligible) == nrow(feature_registry),
          "Complete ledger and feature registry disagree")

annotation_path <- ml_resolve(contract$gene_annotation)
dge_path <- ml_resolve(contract$dge_rds)
signature_path <- file.path(
  ml_source_root(contract), "reproduction", "signature_genes.tsv"
)
ml_assert(file.exists(signature_path), "Missing released signature mapping")
dge <- readRDS(dge_path)
counts <- if (inherits(dge, "DGEList")) dge$counts else dge$counts
ml_assert(!is.null(counts), "Frozen DGE object has no count matrix")
resource_ids <- unique(ml_base_gene_id(rownames(counts)))
ml_assert(length(resource_ids) == contract$gene_family_size,
          "Expression universe is not the frozen 23,370-gene family")
annotation <- fread(annotation_path, select = c("gene_id", "gene_name"))
annotation[, `:=`(
  gene_id_base = ml_base_gene_id(gene_id),
  gene_symbol = toupper(trimws(gene_name))
)]
annotation <- unique(annotation[
  gene_id_base %in% resource_ids & !is.na(gene_symbol) & nzchar(gene_symbol),
  .(gene_id_base, gene_symbol)
])
ml_assert(uniqueN(annotation$gene_id_base) == length(resource_ids),
          "GENCODE v49 annotation does not cover the expression universe")
signature <- fread(signature_path, na.strings = c("", "NA"))
signature_ids <- unique(signature[in_resource == TRUE, gene_id_base])
signature_symbols <- unique(toupper(signature[in_resource == TRUE, gene_symbol]))
ml_assert(length(signature_ids) == contract$signature_observed_size,
          "Observed 145-signature membership drift")

pathway_membership <- rbindlist(pathway_memberships, use.names = TRUE, fill = TRUE)
pathway_membership <- pathway_membership[
  !source_gene %in% signature_symbols & !source_gene %in% signature_ids
]
pathway_membership[, source_weight := 1 / uniqueN(source_gene), by = feature_id]
all_membership <- rbindlist(
  list(pathway_membership, hotspot_memberships), use.names = TRUE, fill = TRUE
)
all_membership <- all_membership[
  !source_gene %in% signature_symbols & !source_gene %in% signature_ids
]
all_membership[, source_gene := toupper(trimws(source_gene))]

symbol_map <- annotation[, .(gene_id_base, source_gene = gene_symbol)]
symbol_membership <- merge(
  all_membership[!grepl("^ENSG[0-9]+", source_gene)], symbol_map,
  by = "source_gene", all.x = TRUE, allow.cartesian = TRUE, sort = FALSE
)
direct_membership <- copy(all_membership[grepl("^ENSG[0-9]+", source_gene)])
direct_membership[, gene_id_base := ml_base_gene_id(source_gene)]
direct_membership[!gene_id_base %in% resource_ids, gene_id_base := NA_character_]
mapped_membership <- rbindlist(
  list(symbol_membership, direct_membership), use.names = TRUE, fill = TRUE
)
mapping_audit <- mapped_membership[, .(
  n_gene_ids = uniqueN(gene_id_base[!is.na(gene_id_base)]),
  mapped = any(!is.na(gene_id_base))
), by = .(feature_id, source_gene, collection)]
mapped_membership <- unique(mapped_membership[
  !is.na(gene_id_base) & !gene_id_base %in% signature_ids,
  .(feature_id, source_gene, collection, source_weight, gene_id_base)
])
mapped_membership[, mapping_multiplicity := uniqueN(gene_id_base),
                  by = .(feature_id, source_gene)]
mapped_membership[, mapped_weight := source_weight / mapping_multiplicity]
mapped_membership[, mapped_weight := mapped_weight / sum(mapped_weight),
                  by = feature_id]
mapped_membership <- mapped_membership[is.finite(mapped_weight) & mapped_weight > 0]
ml_assert(uniqueN(mapped_membership$feature_id) == nrow(feature_registry),
          "At least one geometry-eligible feature has no mapped gene membership")
ml_assert(all(abs(mapped_membership[, sum(mapped_weight), by = feature_id]$V1 - 1) < 1e-10),
          "Mapped feature weights do not sum to one")

feature_levels <- feature_registry$feature_id
gene_levels <- sort(unique(mapped_membership$gene_id_base))
mapped_membership[, feature_index := match(feature_id, feature_levels)]
mapped_membership[, gene_index := match(gene_id_base, gene_levels)]
incidence <- sparseMatrix(
  i = mapped_membership$feature_index,
  j = mapped_membership$gene_index,
  x = mapped_membership$mapped_weight,
  dims = c(length(feature_levels), length(gene_levels)),
  dimnames = list(feature_levels, gene_levels)
)
gene_document_frequency <- Matrix::colSums(incidence != 0)
gene_idf <- log((1 + nrow(incidence)) / (1 + gene_document_frequency)) + 1
weighted_incidence <- incidence %*% Diagonal(x = as.numeric(gene_idf))
row_norm <- sqrt(Matrix::rowSums(weighted_incidence ^ 2))
ml_assert(all(is.finite(row_norm) & row_norm > 0), "Zero-norm geometry row")
weighted_incidence <- Diagonal(x = 1 / row_norm) %*% weighted_incidence

n_svd <- 100L
set.seed(contract$seed)
svd_fit <- irlba(weighted_incidence, nv = n_svd, nu = n_svd,
                 maxit = 1000L, tol = 1e-5)
latent <- sweep(svd_fit$u, 2L, svd_fit$d, "*")
colnames(latent) <- paste0("svd_", seq_len(ncol(latent)))
rownames(latent) <- feature_levels
latent_norm <- sqrt(rowSums(latent ^ 2))
latent_cosine <- latent / latent_norm

k_neighbors <- 30L
maximum_neighbors <- 40L
similarity <- tcrossprod(latent_cosine)
diag(similarity) <- -Inf
neighbors <- matrix(NA_integer_, nrow(similarity), maximum_neighbors)
neighbor_similarity <- matrix(NA_real_, nrow(similarity), maximum_neighbors)
for (index in seq_len(nrow(similarity))) {
  chosen <- order(similarity[index, ], decreasing = TRUE)[seq_len(maximum_neighbors)]
  neighbors[index, ] <- chosen
  neighbor_similarity[index, ] <- similarity[index, chosen]
}

build_mutual_graph <- function(k_value) {
  directed <- data.table(
    from_index = rep(seq_len(nrow(neighbors)), each = k_value),
    to_index = as.vector(t(neighbors[, seq_len(k_value), drop = FALSE])),
    cosine_similarity = as.vector(t(
      neighbor_similarity[, seq_len(k_value), drop = FALSE]
    ))
  )
  directed[, `:=`(
    node_a = pmin(from_index, to_index),
    node_b = pmax(from_index, to_index)
  )]
  edges <- directed[, .(
    directions = .N,
    cosine_similarity = mean(cosine_similarity)
  ), by = .(node_a, node_b)][directions == 2L]
  edges[, edge_type := paste0("mutual_", k_value, "nn")]
  connected <- unique(c(edges$node_a, edges$node_b))
  isolated <- setdiff(seq_len(nrow(feature_registry)), connected)
  if (length(isolated)) {
    bridges <- data.table(
      node_a = pmin(isolated, neighbors[isolated, 1L]),
      node_b = pmax(isolated, neighbors[isolated, 1L]),
      directions = 1L,
      cosine_similarity = neighbor_similarity[isolated, 1L],
      edge_type = "isolated_node_strongest_neighbor_bridge"
    )
    edges <- unique(rbindlist(list(edges, bridges), use.names = TRUE, fill = TRUE),
                    by = c("node_a", "node_b"))
  }
  edges[, `:=`(
    from = feature_levels[node_a],
    to = feature_levels[node_b]
  )]
  graph <- graph_from_data_frame(
    edges[, .(from, to, weight = pmax(cosine_similarity, 1e-8))],
    directed = FALSE,
    vertices = data.frame(name = feature_levels)
  )
  list(edges = edges, graph = graph, n_isolated_bridges = length(isolated))
}

baseline_graph <- build_mutual_graph(k_neighbors)
mutual <- baseline_graph$edges
graph <- baseline_graph$graph
graph_components <- components(graph)
component_summary <- data.table(
  component_raw = seq_along(graph_components$csize),
  n_features = as.integer(graph_components$csize)
)
setorder(component_summary, -n_features, component_raw)
component_summary[, component_id := sprintf("component_%03d", seq_len(.N))]
main_component_raw <- component_summary$component_raw[[1L]]
component_summary[, component_role := fifelse(
  component_raw == main_component_raw,
  "main_system_component", "excluded_microcomponent"
)]
component_summary[component_role == "excluded_microcomponent",
                  microcomponent_id := sprintf("microcomponent_%02d", seq_len(.N))]
main_feature_ids <- V(graph)$name[graph_components$membership == main_component_raw]
main_graph <- induced_subgraph(graph, vids = main_feature_ids)
ml_assert(vcount(main_graph) == 5439L,
          "Expected 5,439 nodes in the main mutual-30NN component")
ml_assert(sum(component_summary$component_role == "excluded_microcomponent") == 3L &&
            sum(component_summary[component_role == "excluded_microcomponent", n_features]) == 6L,
          "Expected three two-node microcomponents outside the system atlas")
set.seed(contract$seed)
leiden <- cluster_leiden(
  main_graph, objective_function = "modularity", weights = E(main_graph)$weight,
  resolution = 1, n_iterations = 10L
)
set.seed(contract$seed)
leiden_repeat <- cluster_leiden(
  main_graph, objective_function = "modularity", weights = E(main_graph)$weight,
  resolution = 1, n_iterations = 10L
)
ml_assert(identical(as.integer(membership(leiden)),
                    as.integer(membership(leiden_repeat))),
          "Seeded Leiden rerun was not deterministic")

topology_sensitivity <- list(data.table(
  sensitivity_id = "baseline_k30_resolution1",
  k = 30L, resolution = 1,
  n_communities = length(leiden), adjusted_rand_vs_baseline = 1,
  outcome_values_read = FALSE
))
for (alternative_k in c(20L, 40L)) {
  alternative_graph <- build_mutual_graph(alternative_k)$graph
  alternative_graph <- induced_subgraph(alternative_graph, vids = main_feature_ids)
  set.seed(contract$seed)
  alternative_leiden <- cluster_leiden(
    alternative_graph, objective_function = "modularity",
    weights = E(alternative_graph)$weight, resolution = 1,
    n_iterations = 10L
  )
  topology_sensitivity[[length(topology_sensitivity) + 1L]] <- data.table(
    sensitivity_id = paste0("k", alternative_k, "_resolution1"),
    k = alternative_k, resolution = 1,
    n_communities = length(alternative_leiden),
    adjusted_rand_vs_baseline = compare(
      membership(leiden), membership(alternative_leiden), method = "adjusted.rand"
    ),
    outcome_values_read = FALSE
  )
}
for (alternative_resolution in c(0.7, 1.3)) {
  set.seed(contract$seed)
  alternative_leiden <- cluster_leiden(
    main_graph, objective_function = "modularity", weights = E(main_graph)$weight,
    resolution = alternative_resolution, n_iterations = 10L
  )
  topology_sensitivity[[length(topology_sensitivity) + 1L]] <- data.table(
    sensitivity_id = paste0("k30_resolution", alternative_resolution),
    k = 30L, resolution = alternative_resolution,
    n_communities = length(alternative_leiden),
    adjusted_rand_vs_baseline = compare(
      membership(leiden), membership(alternative_leiden), method = "adjusted.rand"
    ),
    outcome_values_read = FALSE
  )
}
topology_sensitivity <- rbindlist(topology_sensitivity)

set.seed(contract$seed)
umap_coordinates <- uwot::umap(
  latent, n_neighbors = k_neighbors, n_components = 2L, metric = "cosine",
  n_epochs = 1000L, init = "spectral", spread = 1, min_dist = 0.15,
  n_threads = 1L,
  n_sgd_threads = 1L, verbose = TRUE, seed = contract$seed
)
set.seed(contract$seed)
umap_repeat <- uwot::umap(
  latent, n_neighbors = k_neighbors, n_components = 2L, metric = "cosine",
  n_epochs = 1000L, init = "spectral", spread = 1, min_dist = 0.15,
  n_threads = 1L,
  n_sgd_threads = 1L, verbose = FALSE, seed = contract$seed
)
max_umap_difference <- max(abs(umap_coordinates - umap_repeat))
ml_assert(is.finite(max_umap_difference) && max_umap_difference < 1e-10,
          "Seeded UMAP rerun was not deterministic")

component_by_feature <- setNames(
  graph_components$membership, V(graph)$name
)
leiden_by_feature <- setNames(as.integer(membership(leiden)), names(membership(leiden)))
strength_by_feature <- setNames(
  as.numeric(strength(graph, weights = E(graph)$weight)), V(graph)$name
)
coordinates <- data.table(
  feature_id = feature_levels,
  umap_1 = umap_coordinates[, 1L],
  umap_2 = umap_coordinates[, 2L],
  component_raw = unname(component_by_feature[feature_levels]),
  leiden_raw = unname(leiden_by_feature[feature_levels]),
  graph_strength = unname(strength_by_feature[feature_levels])
)
coordinates <- merge(coordinates, feature_registry, by = "feature_id",
                     all.x = TRUE, sort = FALSE)
coordinates <- merge(
  coordinates,
  component_summary[, .(
    component_raw, component_id, component_role, microcomponent_id
  )],
  by = "component_raw", all.x = TRUE, sort = FALSE
)
raw_summary <- coordinates[component_role == "main_system_component", .(
  n_features = .N,
  centroid_x = median(umap_1), centroid_y = median(umap_2)
), by = leiden_raw]
setorder(raw_summary, -n_features, centroid_x, centroid_y, leiden_raw)
raw_summary[, community_id := sprintf("system_%02d", seq_len(.N))]
coordinates <- merge(
  coordinates, raw_summary[, .(leiden_raw, community_id)],
  by = "leiden_raw", all.x = TRUE, sort = FALSE
)

medoids <- coordinates[component_role == "main_system_component", {
  candidates <- .SD[order(-graph_strength, display_label, feature_id)]
  candidates[1L, .(
    medoid_feature_id = feature_id,
    medoid_feature_label = substr(display_label, 1L, 58L),
    medoid_collection = collection
  )]
}, by = community_id]
community_registry <- coordinates[component_role == "main_system_component", .(
  n_features = .N,
  n_pathways = sum(feature_type == "pathway"),
  n_hotspot_programs = sum(feature_type == "hotspot_program"),
  n_hallmark = sum(collection == "hallmark"),
  n_kegg = sum(collection == "kegg"),
  n_reactome = sum(collection == "reactome"),
  n_go_bp = sum(collection == "go_bp"),
  n_go_mf = sum(collection == "go_mf"),
  n_go_cc = sum(collection == "go_cc"),
  centroid_x = median(umap_1), centroid_y = median(umap_2)
), by = community_id]
community_registry <- merge(community_registry, medoids, by = "community_id",
                            all.x = TRUE, sort = FALSE)
hotspot_summary <- coordinates[
  component_role == "main_system_component" & feature_type == "hotspot_program"
][
  order(community_id, -graph_strength, display_label),
  .(
    representative_hotspot_programs = paste(head(display_label, 5L), collapse = ";"),
    representative_hotspot_cell_types = paste(
      sort(unique(head(cell_type, 5L))), collapse = ";"
    )
  ),
  by = community_id
]
community_registry <- merge(community_registry, hotspot_summary, by = "community_id",
                            all.x = TRUE, sort = FALSE)
community_registry[is.na(representative_hotspot_programs),
                   representative_hotspot_programs := ""]
community_registry[is.na(representative_hotspot_cell_types),
                   representative_hotspot_cell_types := ""]
community_registry[, system_display := community_id]
community_registry[, community_order := as.integer(sub("system_", "", community_id))]
setorder(community_registry, community_order)
community_registry[, community_order := NULL]
coordinates <- merge(
  coordinates,
  community_registry[, .(community_id, system_display)],
  by = "community_id", all.x = TRUE, sort = FALSE
)
coordinates[, feature_order := match(feature_id, feature_levels)]
setorder(coordinates, feature_order)
coordinates[, feature_order := NULL]

community_representatives <- coordinates[component_role == "main_system_component"][
  order(community_id, -graph_strength, display_label, feature_id),
  head(.SD, 10L),
  by = community_id,
  .SDcols = c(
    "feature_id", "collection", "feature_type", "display_label",
    "cell_type", "module", "graph_strength"
  )
]
community_representatives[, representative_rank := seq_len(.N), by = community_id]
setcolorder(
  community_representatives,
  c("community_id", "representative_rank",
    setdiff(names(community_representatives),
            c("community_id", "representative_rank")))
)

latent_table <- as.data.table(latent)
latent_table[, feature_id := feature_levels]
setcolorder(latent_table, c("feature_id", setdiff(names(latent_table), "feature_id")))
mutual_edges <- mutual[, .(
  from, to, cosine_similarity, edge_type
)]
mapped_membership <- merge(
  mapped_membership,
  data.table(gene_id_base = gene_levels, idf = as.numeric(gene_idf)),
  by = "gene_id_base", all.x = TRUE, sort = FALSE
)
mapped_membership[, idf_weighted_membership := mapped_weight * idf]

geometry_parameters <- data.table(
  parameter = c(
    "seed", "feature_count", "complete_feature_family", "gene_count",
    "published_signature_size", "observed_signature_genes_excluded",
    "svd_dimensions", "similarity_metric",
    "knn_k", "graph_rule", "leiden_objective", "leiden_resolution",
    "leiden_iterations", "umap_neighbors", "umap_metric", "umap_min_dist",
    "umap_spread", "umap_epochs", "umap_threads",
    "main_component_nodes", "excluded_microcomponent_nodes",
    "max_seeded_umap_rerun_difference"
  ),
  value = as.character(c(
    contract$seed, nrow(feature_registry), nrow(complete_ledger), length(gene_levels),
    145L, length(signature_ids), n_svd, "cosine", k_neighbors,
    "mutual_30nn_plus_isolate_bridge_main_component_only", "modularity", 1, 10,
    k_neighbors, "cosine", 0.15, 1, 1000, 1L, vcount(main_graph),
    nrow(feature_registry) - vcount(main_graph), max_umap_difference
  ))
)

ml_write_tsv_once(complete_ledger,
                  file.path(source_dir, "complete_feature_testability_ledger.tsv.gz"))
ml_write_tsv_once(feature_registry,
                  file.path(source_dir, "geometry_feature_registry.tsv"))
ml_write_tsv_once(mapping_audit,
                  file.path(source_dir, "gene_mapping_audit.tsv.gz"))
ml_write_tsv_once(mapped_membership,
                  file.path(source_dir, "outcome_blind_gene_membership.tsv.gz"))
ml_write_tsv_once(coordinates,
                  file.path(geometry_dir, "latent_coordinates.tsv.gz"))
ml_write_tsv_once(latent_table,
                  file.path(geometry_dir, "svd100_coordinates.tsv.gz"))
ml_write_tsv_once(mutual_edges,
                  file.path(geometry_dir, "mutual_30nn_edges.tsv.gz"))
ml_write_tsv_once(community_registry,
                  file.path(geometry_dir, "community_registry.tsv"))
ml_write_tsv_once(community_representatives,
                  file.path(geometry_dir, "community_top10_representatives.tsv.gz"))
ml_write_tsv_once(component_summary,
                  file.path(geometry_dir, "graph_component_registry.tsv"))
ml_write_tsv_once(geometry_parameters,
                  file.path(geometry_dir, "geometry_parameters.tsv"))
ml_write_tsv_once(topology_sensitivity,
                  file.path(geometry_dir, "topology_sensitivity.tsv"))

input_rows[[length(input_rows) + 1L]] <- data.table(
  path = normalizePath(c(annotation_path, dge_path, signature_path), mustWork = TRUE),
  input_role = c("gencode_v49_mapping", "expression_universe_only", "signature_exclusion"),
  collection = "shared"
)
input_manifest <- rbindlist(input_rows, use.names = TRUE, fill = TRUE)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
forbidden_pattern <- "stage_vs_continuum|meta_analysis|continuum_models|fibrosis_models|participant_scores|donor_scores"
ml_assert(!any(grepl(forbidden_pattern, input_manifest$path)),
          "Outcome file entered outcome-blind geometry manifest")
input_manifest[, outcome_values_read := FALSE]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "geometry_input_manifest.tsv"))

input_manifest_is_outcome_blind <- !any(grepl(forbidden_pattern, input_manifest$path))
testability_selection_matches <-
  sum(complete_ledger$geometry_eligible) == nrow(feature_registry)
signature_exclusion_passes <-
  !any(mapped_membership$gene_id_base %in% signature_ids)
leiden_rerun_identical <- identical(
  as.integer(membership(leiden)), as.integer(membership(leiden_repeat))
)
complete_ledger_passes <-
  nrow(complete_ledger) == sum(expected_sizes) + contract$program_family_size
component_exclusion_passes <-
  sum(coordinates$component_role == "main_system_component") == 5439L &&
  sum(coordinates$component_role == "excluded_microcomponent") == 6L &&
  nrow(community_registry) == 43L
outcome_blind_audit <- data.table(
  check = c(
    "geometry_reads_no_stage_or_continuum_outcomes",
    "node_selection_uses_testability_only",
    "signature_genes_removed_before_similarity",
    "seeded_leiden_membership_identical",
    "seeded_umap_coordinates_identical",
    "complete_family_ledger_preserved",
    "microcomponents_excluded_from_system_summaries"
  ),
  passed = c(
    input_manifest_is_outcome_blind, testability_selection_matches,
    signature_exclusion_passes, leiden_rerun_identical,
    max_umap_difference < 1e-10, complete_ledger_passes,
    component_exclusion_passes
  ),
  detail = c(
    "Input manifest contains memberships, names, testability, annotation, expression-universe IDs, and signature mapping only",
    "Geometry eligibility requires coverage/testability in both evaluation cohorts; no effect estimate is read",
    paste0(length(signature_ids), " observed signature genes removed by GENCODE v49 ID and symbol"),
    "Leiden repeated with seed 20260817",
    format(max_umap_difference, scientific = TRUE),
    paste0(nrow(complete_ledger), " total features; ", nrow(feature_registry), " geometry eligible"),
    "Three disconnected two-node microcomponents remain explicit but are not named or summarized as molecular systems"
  )
)
ml_assert(all(outcome_blind_audit$passed), "Outcome-blind geometry audit failed")
ml_write_tsv_once(outcome_blind_audit,
                  file.path(provenance_dir, "outcome_blind_audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "geometry_sessionInfo.txt"))

geometry_outputs <- c(
  file.path(geometry_dir, "latent_coordinates.tsv.gz"),
  file.path(geometry_dir, "svd100_coordinates.tsv.gz"),
  file.path(geometry_dir, "mutual_30nn_edges.tsv.gz"),
  file.path(geometry_dir, "community_registry.tsv"),
  file.path(geometry_dir, "community_top10_representatives.tsv.gz"),
  file.path(geometry_dir, "graph_component_registry.tsv"),
  file.path(geometry_dir, "geometry_parameters.tsv"),
  file.path(geometry_dir, "topology_sensitivity.tsv"),
  file.path(source_dir, "complete_feature_testability_ledger.tsv.gz"),
  file.path(source_dir, "geometry_feature_registry.tsv"),
  file.path(source_dir, "outcome_blind_gene_membership.tsv.gz")
)
output_checksums <- data.table(
  path = normalizePath(geometry_outputs, mustWork = TRUE),
  sha256 = vapply(geometry_outputs, ml_sha256, character(1))
)
ml_write_tsv_once(output_checksums,
                  file.path(provenance_dir, "geometry_output_checksums.tsv"))
writeLines("geometry_frozen", file.path(geometry_dir, "GEOMETRY_FROZEN.ok"))
