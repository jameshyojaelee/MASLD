import type {
  CommunityData,
  CommunityEntry,
  GeneGraph,
  LayerMetadata,
  NetworkSearchEntry,
  NetworkNode,
  EdgeLayer,
  EdgeType,
  CommunityInfo,
  CommunityTransition,
  NetworkEdge,
  GeneGraphV2,
  LayerMetadataV2,
  NetworkSearchEntryV2,
  StageNeighborhood,
} from "./network-types";
import { dataUrl } from "./data-base";
import { queryParquet, sqlString } from "./duck";
import { layerColor } from "./palette";

const BASE = dataUrl("network");
const BASE_V2 = dataUrl("network/portal_export_v2");

// ---------------------------------------------------------------------------
// Singleton caches — prevent re-fetching across component renders
// ---------------------------------------------------------------------------

let searchIndexCache: NetworkSearchEntry[] | null = null;
let searchIndexPromise: Promise<NetworkSearchEntry[]> | null = null;

let communityCache: CommunityData | null = null;
let communityPromise: Promise<CommunityData> | null = null;

let layerMetaCache: LayerMetadata[] | null = null;
let layerMetaPromise: Promise<LayerMetadata[]> | null = null;

let globalLayoutCache: NetworkNode[] | null = null;
let globalLayoutPromise: Promise<NetworkNode[]> | null = null;

// ---------------------------------------------------------------------------
// Search index (lightweight — symbol + degree + community)
// ---------------------------------------------------------------------------

export async function loadSearchIndex(): Promise<NetworkSearchEntry[]> {
  if (searchIndexCache) return searchIndexCache;
  if (searchIndexPromise) return searchIndexPromise;

  // The legacy `network/search_index.json` is retired; derive the v1-shaped
  // autocomplete entries from the shipped portal_export_v2 search index.
  searchIndexPromise = loadSearchIndexV2()
    .then((idx) =>
      idx.map((e) => ({
        symbol: e.s,
        ensembl_id: "",
        degree: e.d ?? 0,
        community_macro: e.c_macro_F01 ?? e.c_macro_F34 ?? 0,
        community_label: null,
      }))
    )
    .then((data) => {
      searchIndexCache = data;
      return data;
    });

  return searchIndexPromise;
}

// ---------------------------------------------------------------------------
// Community structure
// ---------------------------------------------------------------------------

export async function loadCommunities(): Promise<CommunityData> {
  if (communityCache) return communityCache;
  if (communityPromise) return communityPromise;

  // The legacy `network/communities.json` is retired. Derive macro entries from
  // the shipped portal_export_v2 F0-F1 community summary — used only for the
  // sidebar's community-name lookup; node COLORING reads `node.community_macro`.
  communityPromise = loadCommunitiesF01()
    .then((p) => {
      const macro: CommunityEntry[] = p.communities.map((c) => ({
        id: c.id,
        label: c.label ?? `Community ${c.id}`,
        size: c.n_genes ?? 0,
        top_genes: [],
        top_pathways: c.top_hallmark ? [c.top_hallmark] : [],
        color: "",
      }));
      const data: CommunityData = { macro, meso: [] };
      communityCache = data;
      return data;
    })
    .catch(() => {
      const data: CommunityData = { macro: [], meso: [] };
      communityCache = data;
      return data;
    });

  return communityPromise;
}

// ---------------------------------------------------------------------------
// Per-gene neighborhood graph — see v2 `loadGeneGraph` below. v1 available as
// `v1_loadGeneGraph`.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Global FA2 layout for community map mode
// ---------------------------------------------------------------------------

export async function loadGlobalLayout(): Promise<NetworkNode[]> {
  if (globalLayoutCache) return globalLayoutCache;
  if (globalLayoutPromise) return globalLayoutPromise;

  globalLayoutPromise = fetch(`${BASE}/global_layout.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load global layout: ${r.status}`);
      return r.json() as Promise<NetworkNode[]>;
    })
    .then((data) => {
      globalLayoutCache = data;
      return data;
    });

  return globalLayoutPromise;
}

// ---------------------------------------------------------------------------
// Layer metadata (colors, descriptions, edge counts)
// ---------------------------------------------------------------------------

/** Default layer metadata — used as fallback if layers.json is unavailable.
 *
 * Labels aligned to v4 semantics: `ppi` is the STRING generic prior; the
 * D-prefixed layers are the MASLD-specific dynamic evidence edges.
 * Layers flagged `legacy: true` are retained only for v1 backcompat — the
 * v4 UI filters them out of display.
 */
export const DEFAULT_LAYERS: LayerMetadata[] = [
  { id: "ppi", label: "STRING (generic prior)", color: "#e74c3c", description: "STRING v12 \u2265700 protein-protein backbone (generic knowledge prior, not MASLD-specific)", edge_count: 115020 },
  { id: "coexpr", label: "D-F2 stage-dynamic coexpression", color: "#3498db", description: "Stage-stratified coexpression edges with |\u0394r| FDR<0.05 across F0-F1 / F2 / F3-F4", edge_count: 3260461 },
  { id: "regulon", label: "Regulon", color: "#2ecc71", description: "SCENIC+ TF-target regulons", edge_count: 0, legacy: true },
  { id: "lr", label: "D-LR (F2-differential ligand-receptor)", color: "#f39c12", description: "LIANA F2-differentially rewired ligand-receptor pairs", edge_count: 9 },
  { id: "genetic", label: "D-COLOC (GWAS colocalization)", color: "#9b59b6", description: "SuSiE-COLOC locus pairs with PP.H4\u22650.5 across 17 MASLD GWAS", edge_count: 111 },
  { id: "pathway", label: "Pathway", color: "#1abc9c", description: "Shared pathway membership (MSigDB)", edge_count: 0, legacy: true },
  { id: "spatial", label: "Spatial", color: "#e67e22", description: "Spatial co-localization (Visium)", edge_count: 0, legacy: true },
  { id: "cosmos", label: "COSMOS", color: "#34495e", description: "COSMOS mechanistic paths", edge_count: 0, legacy: true },
  { id: "cerna", label: "D-ceRNA", color: "#e91e63", description: "MASLD-expressed lncRNA ceRNA triplets", edge_count: 0 },
  { id: "xspecies", label: "D-XSpecies", color: "#607d8b", description: "Cross-species ortholog early-to-late transition pairs (FPC model)", edge_count: 0 },
];

export async function loadLayerMetadata(): Promise<LayerMetadata[]> {
  if (layerMetaCache) return layerMetaCache;
  if (layerMetaPromise) return layerMetaPromise;

  layerMetaPromise = fetch(`${BASE}/layers.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load layer metadata: ${r.status}`);
      return r.json() as Promise<LayerMetadata[]>;
    })
    .then((data) => {
      layerMetaCache = data;
      return data;
    })
    .catch(() => {
      // Fallback to hardcoded defaults if file not yet deployed
      layerMetaCache = DEFAULT_LAYERS;
      return DEFAULT_LAYERS;
    });

  return layerMetaPromise;
}

// ---------------------------------------------------------------------------
// Utility: get layer color by id
// ---------------------------------------------------------------------------

// Layer→color is owned by the palette authority (`LAYER_COLORS`/`layerColor`),
// so the canvas, chart legends, and controls agree and the gray fallback is the
// canonical control gray.
export function getLayerColor(layer: EdgeLayer): string {
  return layerColor(layer);
}

// ===========================================================================
// v2 loaders (Architecture C) — portal_export_v2
// ===========================================================================
//
// Backward-compat: legacy loaders above are unchanged. Components that need
// v1 payloads can import them under the `v1_` alias re-exports at the bottom.

interface CommunitiesV2Payload {
  side: "F01" | "F34";
  n_communities: number;
  communities: CommunityInfo[];
}

interface TransitionsV2Payload {
  n_transitions: number;
  transitions: CommunityTransition[];
}

let communitiesF01Cache: CommunitiesV2Payload | null = null;
let communitiesF01Promise: Promise<CommunitiesV2Payload> | null = null;
let communitiesF34Cache: CommunitiesV2Payload | null = null;
let communitiesF34Promise: Promise<CommunitiesV2Payload> | null = null;
let transitionsCache: TransitionsV2Payload | null = null;
let transitionsPromise: Promise<TransitionsV2Payload> | null = null;

let layerMetaV2Cache: LayerMetadataV2[] | null = null;
let layerMetaV2Promise: Promise<LayerMetadataV2[]> | null = null;

let searchIndexV2Cache: NetworkSearchEntryV2[] | null = null;
let searchIndexV2Promise: Promise<NetworkSearchEntryV2[]> | null = null;

const geneGraphV2Cache = new Map<string, GeneGraphV2>();

export async function loadCommunitiesF01(): Promise<CommunitiesV2Payload> {
  if (communitiesF01Cache) return communitiesF01Cache;
  if (communitiesF01Promise) return communitiesF01Promise;
  communitiesF01Promise = fetch(`${BASE_V2}/communities_F01.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load communities_F01: ${r.status}`);
      return r.json() as Promise<CommunitiesV2Payload>;
    })
    .then((data) => {
      communitiesF01Cache = data;
      return data;
    });
  return communitiesF01Promise;
}

export async function loadCommunitiesF34(): Promise<CommunitiesV2Payload> {
  if (communitiesF34Cache) return communitiesF34Cache;
  if (communitiesF34Promise) return communitiesF34Promise;
  communitiesF34Promise = fetch(`${BASE_V2}/communities_F34.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load communities_F34: ${r.status}`);
      return r.json() as Promise<CommunitiesV2Payload>;
    })
    .then((data) => {
      communitiesF34Cache = data;
      return data;
    });
  return communitiesF34Promise;
}

export async function loadCommunityTransitions(): Promise<TransitionsV2Payload> {
  if (transitionsCache) return transitionsCache;
  if (transitionsPromise) return transitionsPromise;
  transitionsPromise = fetch(`${BASE_V2}/community_transitions.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load community_transitions: ${r.status}`);
      return r.json() as Promise<TransitionsV2Payload>;
    })
    .then((data) => {
      transitionsCache = data;
      return data;
    });
  return transitionsPromise;
}

export async function loadLayerMetadataV2(): Promise<LayerMetadataV2[]> {
  if (layerMetaV2Cache) return layerMetaV2Cache;
  if (layerMetaV2Promise) return layerMetaV2Promise;
  layerMetaV2Promise = fetch(`${BASE_V2}/layer_metadata.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load layer_metadata v2: ${r.status}`);
      return r.json() as Promise<{ layers: LayerMetadataV2[] }>;
    })
    .then((data) => {
      layerMetaV2Cache = data.layers;
      return data.layers;
    });
  return layerMetaV2Promise;
}

export async function loadSearchIndexV2(): Promise<NetworkSearchEntryV2[]> {
  if (searchIndexV2Cache) return searchIndexV2Cache;
  if (searchIndexV2Promise) return searchIndexV2Promise;
  searchIndexV2Promise = fetch(`${BASE_V2}/search_index.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load search index v2: ${r.status}`);
      return r.json() as Promise<{ genes: NetworkSearchEntryV2[] }>;
    })
    .then((data) => {
      searchIndexV2Cache = data.genes;
      return data.genes;
    });
  return searchIndexV2Promise;
}

/**
 * v2 per-gene annotation loader.
 *
 * The 16,601 `portal_export_v2/gene_graphs/*.json` files are retired (dropped
 * from the HF dataset). The force graph now builds its subgraph from
 * `network_edges.parquet` via `loadGeneGraph` below; this loader survives only
 * to feed the sidebar's community-membership panel, which reads `.community`.
 * We reconstruct that from the shipped v2 search index (per-gene community ids)
 * plus the F0-F1 / F3-F4 community summaries (for the hallmark labels). The
 * heavy `neighbors` / `edge_counts` payload is left empty (unused by the page).
 */
export async function loadGeneGraphV2(symbol: string): Promise<GeneGraphV2> {
  const key = symbol.toUpperCase();
  const cached = geneGraphV2Cache.get(key);
  if (cached) return cached;

  const [idx, cF01, cF34] = await Promise.all([
    loadSearchIndexV2().catch(() => [] as NetworkSearchEntryV2[]),
    loadCommunitiesF01().then((p) => p.communities).catch(() => [] as CommunityInfo[]),
    loadCommunitiesF34().then((p) => p.communities).catch(() => [] as CommunityInfo[]),
  ]);
  const entry = idx.find((e) => e.s.toUpperCase() === key);
  const hallmarkOf = (list: CommunityInfo[], id: number | undefined): string | null =>
    id == null ? null : list.find((c) => c.id === id)?.top_hallmark ?? null;

  const emptyNeighbors: StageNeighborhood = {
    neighbors_string: [],
    neighbors_d_f2_emerging: [],
    neighbors_d_f2_dissolving: [],
    neighbors_d_coloc: [],
    neighbors_d_lr: [],
    top_d_cerna: [],
    top_d_xs: [],
    contested: [],
    druggable_pairs: [],
  };

  const graph: GeneGraphV2 = {
    gene: key,
    ensembl_id: "",
    biotype: null,
    attributes: {
      is_deg: entry?.is_deg,
      dgidb_druggable: entry?.is_drug,
    },
    community: {
      f01:
        entry?.c_macro_F01 != null
          ? {
              macro_id: entry.c_macro_F01,
              community_id: entry.c_macro_F01,
              top_hallmark: hallmarkOf(cF01, entry.c_macro_F01),
            }
          : null,
      f34:
        entry?.c_macro_F34 != null
          ? {
              macro_id: entry.c_macro_F34,
              community_id: entry.c_macro_F34,
              top_hallmark: hallmarkOf(cF34, entry.c_macro_F34),
            }
          : null,
    },
    neighbors: emptyNeighbors,
    edge_counts: {} as Record<EdgeType, number>,
    f_stage_trajectory: { r_F01: null, r_F2: null, r_F34: null, emergence_stage: null },
  };
  geneGraphV2Cache.set(key, graph);
  return graph;
}

/**
 * Legacy v1 per-gene graph loader. Kept for the force-directed rendering path
 * which still consumes the v1 `GeneGraph { nodes, links }` shape.
 */
const v1_geneGraphCache = new Map<string, GeneGraph>();

// Per-gene incident-edge cap. react-force-graph-2d stutters past a few hundred
// edges/gene; the source parquet already applied a per-node top-K, but a hub can
// still accumulate >1k incident edges (it is in many partners' top-K), so cap by
// edge strength here.
const NEIGHBOR_CAP = 300;

type V4Type = NonNullable<NetworkEdge["v4"]>["type"];

/** Raw edge row projected out of `network_edges.parquet`. */
interface EdgeRow {
  gene_a: string;
  gene_b: string;
  type: string | null;
  in_string_ge700: boolean | null;
  is_contested: boolean | null;
  string_score: number | null;
  delta_max: number | null;
  emergence_stage: string | null;
  loco_replication_fraction: number | null;
  pp4_min: number | null;
  gwas_list: string | null;
  score_diff_lr: number | null;
  cell_type_pairs: string | null;
  coloc_linked: boolean | null;
  druggable_pair: boolean | null;
}

/** Raw node row projected out of `network_nodes.parquet`. */
interface NodeRow {
  symbol: string;
  ensembl_id: string | null;
  gene_biotype: string | null;
  bulk_logFC: number | null;
  bulk_padj: number | null;
  is_deg: boolean | null;
  is_conserved_core: boolean | null;
  sex_class: string | null;
  dgidb_druggable: boolean | null;
}

const numOrNull = (x: unknown): number | null =>
  typeof x === "number" && !Number.isNaN(x) ? x : null;
const clamp01 = (x: number): number => Math.max(0, Math.min(1, x));

/**
 * Reconstruct the UI's (edge layer, v4 type, paint posterior) from the parquet's
 * collapsed `type` channel plus the richer signal columns. The parquet stores one
 * primary channel per edge (coexpr/coloc/lr/regulon/cerna/xspecies); within
 * `coexpr` we split STRING backbone (S -> layer `ppi`) from stage-dynamic
 * coexpression (D-F2 -> layer `coexpr`) using `delta_max` / `in_string_ge700`.
 */
function classifyEdge(row: EdgeRow): {
  layer: EdgeLayer;
  v4type: V4Type | undefined;
  posterior: number;
} {
  switch (row.type) {
    case "coloc":
      return { layer: "genetic", v4type: "D-COLOC", posterior: numOrNull(row.pp4_min) ?? 0.5 };
    case "lr":
      return { layer: "lr", v4type: "D-LR", posterior: 0.5 };
    case "regulon":
      return { layer: "regulon", v4type: undefined, posterior: 0.5 };
    case "cerna":
      return { layer: "cerna", v4type: "D-ceRNA", posterior: 0.5 };
    case "xspecies":
      return { layer: "xspecies", v4type: "D-XS", posterior: 0.5 };
    default:
      // "coexpr" (or unknown): stage-dynamic edge wins over the STRING backbone.
      if (numOrNull(row.delta_max) !== null) {
        return {
          layer: "coexpr",
          v4type: "D-F2",
          posterior: numOrNull(row.loco_replication_fraction) ?? 0.5,
        };
      }
      if (row.in_string_ge700) {
        return { layer: "ppi", v4type: "S", posterior: numOrNull(row.string_score) ?? 0.5 };
      }
      return { layer: "coexpr", v4type: undefined, posterior: 0.5 };
  }
}

// symbol -> {community, degree} derived from the shipped v2 search index. The
// network parquet's `community_macro_*` cols are null, so node COLORING and hub
// sizing are sourced from portal_export_v2 (data contract Section 6.2 note).
let nodeIndexCache: Map<string, { community: number; degree: number }> | null = null;
let nodeIndexPromise: Promise<Map<string, { community: number; degree: number }>> | null = null;

async function getNodeIndexMap(): Promise<Map<string, { community: number; degree: number }>> {
  if (nodeIndexCache) return nodeIndexCache;
  if (nodeIndexPromise) return nodeIndexPromise;
  nodeIndexPromise = loadSearchIndexV2()
    .then((idx) => {
      const m = new Map<string, { community: number; degree: number }>();
      for (const e of idx) {
        m.set(e.s.toUpperCase(), {
          community: e.c_macro_F01 ?? e.c_macro_F34 ?? 0,
          degree: e.d ?? 0,
        });
      }
      nodeIndexCache = m;
      return m;
    })
    .catch(() => {
      const m = new Map<string, { community: number; degree: number }>();
      nodeIndexCache = m;
      return m;
    });
  return nodeIndexPromise;
}

/**
 * Build a gene's 1-hop neighborhood subgraph from `network_edges.parquet`
 * (`WHERE gene_a = ? OR gene_b = ?`) + node attributes from
 * `network_nodes.parquet`, via DuckDB-WASM. Returns the same
 * `GeneGraph { center, nodes, links }` shape the force graph consumed from the
 * retired per-gene JSONs, so `force-graph.tsx` is unchanged.
 */
export async function loadGeneGraph(symbol: string): Promise<GeneGraph> {
  const key = symbol.toUpperCase();
  const cached = v1_geneGraphCache.get(key);
  if (cached) return cached;

  const [edgeRows, nodeIdx] = await Promise.all([
    queryParquet<EdgeRow>(
      "network_edges.parquet",
      (t) =>
        `SELECT gene_a, gene_b, type, in_string_ge700, is_contested, string_score,
                delta_max, emergence_stage, loco_replication_fraction, pp4_min,
                gwas_list, score_diff_lr, cell_type_pairs, coloc_linked, druggable_pair,
                GREATEST(COALESCE(ABS(delta_max), 0), COALESCE(pp4_min, 0),
                         COALESCE(string_score, 0),
                         CASE WHEN coloc_linked THEN 0.8 ELSE 0 END,
                         CASE WHEN is_contested THEN 0.6 ELSE 0 END) AS _strength
         FROM ${t}
         WHERE gene_a = ? OR gene_b = ?
         ORDER BY _strength DESC
         LIMIT ${NEIGHBOR_CAP}`,
      [key, key]
    ),
    getNodeIndexMap(),
  ]);

  const neighborSyms = new Set<string>();
  for (const r of edgeRows) {
    const a = (r.gene_a ?? "").toUpperCase();
    const b = (r.gene_b ?? "").toUpperCase();
    const other = a === key ? b : a;
    if (other && other !== key) neighborSyms.add(other);
  }

  const allSyms = [key, ...neighborSyms];
  const nodeAttrs = allSyms.length
    ? await queryParquet<NodeRow>(
        "network_nodes.parquet",
        (t) =>
          `SELECT symbol, ensembl_id, gene_biotype, bulk_logFC, bulk_padj,
                  is_deg, is_conserved_core, sex_class, dgidb_druggable
           FROM ${t}
           WHERE symbol IN (${allSyms.map(sqlString).join(",")})`
      )
    : [];
  const attrBySym = new Map<string, NodeRow>();
  for (const a of nodeAttrs) attrBySym.set(String(a.symbol).toUpperCase(), a);

  // A gene with no incident edges AND no node row is genuinely absent from the
  // network — throw so the page surfaces "Could not load network for …".
  if (edgeRows.length === 0 && !attrBySym.has(key)) {
    throw new Error(`Gene not in network: ${key}`);
  }

  const mkNode = (sym: string): NetworkNode => {
    const a = attrBySym.get(sym);
    const idx = nodeIdx.get(sym);
    return {
      id: sym,
      symbol: sym,
      ensembl_id: (a?.ensembl_id as string | null) ?? "",
      biotype: (a?.gene_biotype as string | null) ?? null,
      community_macro: idx?.community ?? 0,
      community_meso: 0,
      community_micro: 0,
      community_label: null,
      degree: idx?.degree ?? (sym === key ? edgeRows.length : 1),
      layers_active: 0,
      evidence: {} as NetworkNode["evidence"],
      bulk_logfc: numOrNull(a?.bulk_logFC),
      bulk_padj: numOrNull(a?.bulk_padj),
      sex_class: (a?.sex_class as string | null) ?? null,
      progression_class: null,
      is_deg: !!a?.is_deg,
      is_conserved_core: !!a?.is_conserved_core,
      dgidb_druggable: !!a?.dgidb_druggable,
      fa2_x: null,
      fa2_y: null,
    };
  };

  const nodeMap = new Map<string, NetworkNode>();
  nodeMap.set(key, mkNode(key));
  for (const s of neighborSyms) if (!nodeMap.has(s)) nodeMap.set(s, mkNode(s));

  // Star topology (center -> each partner), matching the retired per-gene export.
  const links: NetworkEdge[] = [];
  for (const r of edgeRows) {
    const a = (r.gene_a ?? "").toUpperCase();
    const b = (r.gene_b ?? "").toUpperCase();
    const partner = a === key ? b : a;
    if (!partner || partner === key || !nodeMap.has(partner)) continue;
    const { layer, v4type, posterior } = classifyEdge(r);
    links.push({
      source: key,
      target: partner,
      layer,
      posterior: clamp01(posterior),
      color: layerColor(layer),
      v4: v4type
        ? {
            type: v4type,
            emergence_stage: r.emergence_stage ?? undefined,
            loco_replication: numOrNull(r.loco_replication_fraction) ?? undefined,
            pp4_min: numOrNull(r.pp4_min) ?? undefined,
            druggable_pair: !!r.druggable_pair,
          }
        : undefined,
    });
  }

  // Fall back to incident-edge count for the center's size when the search index
  // has no global degree for it.
  const center = nodeMap.get(key)!;
  if (!nodeIdx.get(key)) center.degree = links.length;

  const graph: GeneGraph = {
    center: key,
    nodes: Array.from(nodeMap.values()),
    links,
  };
  v1_geneGraphCache.set(key, graph);
  return graph;
}

export const v1_loadGeneGraph = loadGeneGraph;
export const v1_loadSearchIndex = loadSearchIndex;
export const v1_loadCommunities = loadCommunities;
export const v1_loadLayerMetadata = loadLayerMetadata;
export const v1_loadGlobalLayout = loadGlobalLayout;
