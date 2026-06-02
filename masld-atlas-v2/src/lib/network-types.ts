import type { EvidenceStrengths } from "./types";

/** Individual node in the force graph */
export interface NetworkNode {
  id: string;
  symbol: string;
  ensembl_id: string;
  biotype: string | null;
  community_macro: number;
  community_meso: number;
  community_micro: number;
  community_label: string | null;
  degree: number;
  layers_active: number;
  evidence: EvidenceStrengths;
  dream_logfc: number | null;
  dream_padj: number | null;
  sex_class: string | null;
  progression_class: string | null;
  is_deg: boolean;
  is_conserved_core: boolean;
  dgidb_druggable: boolean;
  /** Pre-computed FA2 layout x (community map mode) */
  fa2_x: number | null;
  /** Pre-computed FA2 layout y (community map mode) */
  fa2_y: number | null;
  /** Populated at runtime by d3-force */
  x?: number;
  y?: number;
  vx?: number;
  vy?: number;
}

/** Edge between two nodes */
export interface NetworkEdge {
  source: string;
  target: string;
  layer: EdgeLayer;
  posterior: number;
  /** Hex color derived from layer */
  color: string;
  /**
   * v4 edge metadata passthrough (optional — only populated for adapters that
   * load from `portal_export_v2/gene_graphs/*.json`). Fields mirror
   * `NetworkEdgeV2`, but only the subset the adapter can confidently extract.
   */
  v4?: {
    type: "S" | "D-F2" | "D-COLOC" | "D-LR" | "D-ceRNA" | "D-XS";
    emergence_stage?: string;
    loco_replication?: number;
    pp4_min?: number;
    druggable_pair?: boolean;
    conserved_mouse_a?: boolean;
    conserved_mouse_b?: boolean;
  };
}

/** The 10 edge modality layers */
export type EdgeLayer =
  | "ppi"
  | "coexpr"
  | "regulon"
  | "lr"
  | "genetic"
  | "pathway"
  | "spatial"
  | "cosmos"
  | "cerna"
  | "xspecies";

/** Per-gene subgraph loaded on demand */
export interface GeneGraph {
  center: string;
  nodes: NetworkNode[];
  links: NetworkEdge[];
}

/** Metadata for a single edge layer */
export interface LayerMetadata {
  id: EdgeLayer;
  label: string;
  color: string;
  description: string;
  edge_count: number;
  /** Marks layers that are v1-legacy (retained for backcompat, hidden in v4 UI). */
  legacy?: boolean;
}

/** A macro community entry from communities.json */
export interface CommunityEntry {
  id: number;
  label: string;
  size: number;
  top_genes: string[];
  top_pathways: string[];
  color: string;
}

/** Full community data structure */
export interface CommunityData {
  macro: CommunityEntry[];
  meso: CommunityEntry[];
}

/** Halo data attached to a gene (drugs, variants, pathways) */
export interface GeneHalo {
  drugs: string[];
  gwas_variants: string[];
  pathways: string[];
  regulons: string[];
}

/** Entry in the search index (lightweight, for autocomplete) */
export interface NetworkSearchEntry {
  symbol: string;
  ensembl_id: string;
  degree: number;
  community_macro: number;
  community_label: string | null;
}

/** View modes for the network explorer */
export type NetworkViewMode = "neighborhood" | "community" | "compare";

// ===========================================================================
// v2 (Architecture C) types
// ===========================================================================

/**
 * v2 edge types: STRING baseline (S) plus 5 D-type evidence layers.
 * Replaces the legacy 10-layer `EdgeLayer` union for v2 portals.
 */
export type EdgeType =
  | "S"
  | "D-F2"
  | "D-COLOC"
  | "D-LR"
  | "D-ceRNA"
  | "D-XS";

/**
 * F2-switch emergence classification for an edge or gene.
 * `null` indicates the assignment was not computable (e.g. insufficient stage coverage).
 */
export type EmergenceStage =
  | "F2_emerging"
  | "F2_dissolving"
  | "progressive_up"
  | "progressive_down"
  | "transient_F2"
  | "F34_specific"
  | "invariant"
  | null;

/**
 * Full edge annotation record from `edge_annotation_atlas.parquet`, oriented
 * at a specific query gene (the `partner` is the non-query endpoint).
 * Optional fields are populated only for edges of relevant type.
 */
export interface NetworkEdgeV2 {
  partner: string;
  type: EdgeType;
  /** STRING combined score (v12, 0-1 normalized); present for S-type and contested edges. */
  string_score?: number;
  in_string_ge700?: boolean;
  /** Edge is present in both STRING and at least one D-type layer. */
  is_contested?: boolean;

  // Stage-dynamic co-expression (D-F2)
  r_F01?: number;
  r_F2?: number;
  r_F34?: number;
  delta_max?: number;
  emergence_stage?: EmergenceStage;
  fdr_f2?: number;
  loco_replication_fraction?: number;

  // Genetic colocalization (D-COLOC)
  pp4_min?: number;
  gwas_list?: string;
  locus_id?: string;

  // Ligand-receptor (D-LR)
  score_diff_lr?: number;
  cell_type_pairs?: string;
  emergence_stage_lr?: EmergenceStage;
  celltype_driver_a?: string;
  celltype_driver_b?: string;

  // Provenance cross-links
  coloc_linked_gene_a?: string;
  coloc_linked_gene_b?: string;
  coloc_linked?: boolean;
  druggable_pair?: boolean;
  druggable_a?: boolean;
  druggable_b?: boolean;

  // Endpoint annotations
  sex_class_a?: string;
  sex_class_b?: string;
  conserved_mouse_a?: boolean;
  conserved_mouse_b?: boolean;
  ferroptosis_a?: string;
  ferroptosis_b?: string;
  zonation_a?: string;
  zonation_b?: string;

  provenance_hash?: string;
}

/** Stage-specific neighborhood buckets emitted by the v2 portal export. */
export interface StageNeighborhood {
  neighbors_string: NetworkEdgeV2[];
  neighbors_d_f2_emerging: NetworkEdgeV2[];
  neighbors_d_f2_dissolving: NetworkEdgeV2[];
  neighbors_d_coloc: NetworkEdgeV2[];
  neighbors_d_lr: NetworkEdgeV2[];
  top_d_cerna: NetworkEdgeV2[];
  top_d_xs: NetworkEdgeV2[];
  contested: NetworkEdgeV2[];
  druggable_pairs: NetworkEdgeV2[];
}

/** Gene-level attributes surfaced in the v2 per-gene JSON. */
export interface GeneAttributesV2 {
  is_deg?: boolean;
  dream_logFC?: number;
  dream_padj?: number;
  coloc_susie_best_pp4?: number;
  is_conserved_core?: boolean;
  sex_class?: string;
  dgidb_druggable?: boolean;
  ferroptosis_class?: string | null;
  zonation_class?: string | null;
  attribution_class?: string | null;
}

/** Community membership entry (F01 or F34 side). */
export interface CommunityMembership {
  community_id?: number;
  macro_id?: number;
  meso_id?: number;
  micro_id?: number;
  top_hallmark?: string | null;
  top_hallmark_p?: number | null;
}

/** A single community with its pathway enrichment and dominance metrics. */
export interface CommunityInfo {
  id: number;
  label: string;
  n_genes: number;
  top_hallmark: string | null;
  top_hallmark_p: number | null;
  overlap?: number | null;
  set_size?: number | null;
  dominance_change: number | null;
  best_partner_other_side?: number | null;
  best_shared?: number | null;
}

/** An F0-F1 ↔ F3-F4 community transition entry. */
export interface CommunityTransition {
  f01_id: number;
  f34_id: number;
  n_shared: number;
  n_F01: number;
  n_F34: number;
  jaccard: number;
  F01_class: string | null;
  F34_class: string | null;
  classification: string | null;
}

/** v2 per-gene graph document (from `portal_export_v2/gene_graphs/{SYMBOL}.json`). */
export interface GeneGraphV2 {
  gene: string;
  ensembl_id?: string;
  biotype?: string | null;
  attributes: GeneAttributesV2;
  community: {
    f01: CommunityMembership | null;
    f34: CommunityMembership | null;
  };
  neighbors: StageNeighborhood;
  edge_counts: Record<EdgeType, number>;
  f_stage_trajectory: {
    r_F01: number | null;
    r_F2: number | null;
    r_F34: number | null;
    emergence_stage: EmergenceStage;
  };
}

/** v2 layer metadata entry. */
export interface LayerMetadataV2 {
  type: EdgeType;
  color: string;
  description: string;
  n_edges: number;
  gold_fold_enrichment: number | null;
}

/** v2 compact search-index entry (one per gene). */
export interface NetworkSearchEntryV2 {
  s: string;
  d: number;
  c_macro_F01?: number;
  c_macro_F34?: number;
  is_deg?: boolean;
  is_coloc?: boolean;
  is_drug?: boolean;
}

/** Node coloring scheme */
export type NodeColorBy =
  | "community"
  | "layers"
  | "sex_class"
  | "progression"
  | "logfc"
  | "druggability";
