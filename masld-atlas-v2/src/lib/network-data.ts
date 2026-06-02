import type {
  CommunityData,
  GeneGraph,
  LayerMetadata,
  NetworkSearchEntry,
  NetworkNode,
  EdgeLayer,
  CommunityInfo,
  CommunityTransition,
  NetworkEdge,
  GeneGraphV2,
  LayerMetadataV2,
  NetworkSearchEntryV2,
} from "./network-types";

const BASE = "/data/network";
const BASE_V2 = "/data/network/portal_export_v2";

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

  searchIndexPromise = fetch(`${BASE}/search_index.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load search index: ${r.status}`);
      return r.json() as Promise<NetworkSearchEntry[]>;
    })
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

  communityPromise = fetch(`${BASE}/communities.json`)
    .then((r) => {
      if (!r.ok) throw new Error(`Failed to load communities: ${r.status}`);
      return r.json() as Promise<CommunityData>;
    })
    .then((data) => {
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
  { id: "xspecies", label: "D-XSpecies", color: "#607d8b", description: "Cross-species ortholog F2 switch pairs (FPC model)", edge_count: 0 },
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

const LAYER_COLOR_MAP: Record<EdgeLayer, string> = {
  ppi: "#e74c3c",
  coexpr: "#3498db",
  regulon: "#2ecc71",
  lr: "#f39c12",
  genetic: "#9b59b6",
  pathway: "#1abc9c",
  spatial: "#e67e22",
  cosmos: "#34495e",
  cerna: "#e91e63",
  xspecies: "#607d8b",
};

export function getLayerColor(layer: EdgeLayer): string {
  return LAYER_COLOR_MAP[layer] ?? "#888888";
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
 * v2 per-gene graph loader. Fetches from `portal_export_v2/gene_graphs/`.
 * Throws if the gene graph is unavailable.
 */
export async function loadGeneGraphV2(symbol: string): Promise<GeneGraphV2> {
  const key = symbol.toUpperCase();
  const cached = geneGraphV2Cache.get(key);
  if (cached) return cached;

  const response = await fetch(`${BASE_V2}/gene_graphs/${key}.json`);
  if (!response.ok) {
    throw new Error(`Gene graph v2 not found: ${key} (${response.status})`);
  }
  const data = (await response.json()) as GeneGraphV2;
  geneGraphV2Cache.set(key, data);
  return data;
}

/**
 * Legacy v1 per-gene graph loader. Kept for the force-directed rendering path
 * which still consumes the v1 `GeneGraph { nodes, links }` shape.
 */
const v1_geneGraphCache = new Map<string, GeneGraph>();

// Map v4 edge type -> v1 EdgeLayer so the legacy force-graph renderer can
// consume portal_export_v2/gene_graphs/*.json unchanged.
const V4_TYPE_TO_LAYER: Record<string, EdgeLayer> = {
  S: "ppi",
  "D-F2": "coexpr",
  "D-COLOC": "genetic",
  "D-LR": "lr",
  "D-ceRNA": "cerna",
  "D-XS": "xspecies",
};

export async function loadGeneGraph(symbol: string): Promise<GeneGraph> {
  const key = symbol.toUpperCase();
  const cached = v1_geneGraphCache.get(key);
  if (cached) return cached;

  const response = await fetch(`${BASE_V2}/gene_graphs/${key}.json`);
  if (!response.ok) {
    throw new Error(`Gene graph not found: ${key} (${response.status})`);
  }
  const v2 = (await response.json()) as GeneGraphV2;

  const nodeMap = new Map<string, NetworkNode>();
  const mkNode = (sym: string): NetworkNode => ({
    id: sym,
    symbol: sym,
    ensembl_id: "",
    biotype: null,
    community_macro: 0,
    community_meso: 0,
    community_micro: 0,
    community_label: null,
    degree: 0,
    layers_active: 0,
    evidence: {} as NetworkNode["evidence"],
    dream_logfc: null,
    dream_padj: null,
    sex_class: null,
    progression_class: null,
    is_deg: false,
    is_conserved_core: false,
    dgidb_druggable: false,
    fa2_x: null,
    fa2_y: null,
  });
  nodeMap.set(key, mkNode(key));

  const links: NetworkEdge[] = [];
  const buckets: Array<keyof GeneGraphV2["neighbors"]> = [
    "neighbors_string",
    "neighbors_d_f2_emerging",
    "neighbors_d_f2_dissolving",
    "neighbors_d_coloc",
    "neighbors_d_lr",
    "top_d_cerna",
    "top_d_xs",
  ];
  for (const b of buckets) {
    const arr = (v2.neighbors?.[b] ?? []) as Array<{
      partner?: string;
      type?: string;
      string_score?: number;
      emergence_stage?: string;
      emergence_stage_lr?: string;
      loco_replication_fraction?: number;
      pp4_min?: number;
      druggable_pair?: boolean;
      conserved_mouse_a?: boolean;
      conserved_mouse_b?: boolean;
    }>;
    for (const e of arr) {
      const partner = (e.partner ?? "").toUpperCase();
      if (!partner) continue;
      if (!nodeMap.has(partner)) nodeMap.set(partner, mkNode(partner));
      const t = e.type ?? "";
      const layer = V4_TYPE_TO_LAYER[t] ?? "coexpr";
      const v4Type = (["S", "D-F2", "D-COLOC", "D-LR", "D-ceRNA", "D-XS"] as const).find(
        (x) => x === t
      );
      // Posterior semantics differ per edge type — we stash the type-appropriate
      // strength here for the v1 link-paint gradient, but the edgeVisible() path
      // in force-graph.tsx consults `v4` for the real threshold dispatch.
      let posterior = 0.5;
      if (t === "S" && typeof e.string_score === "number") posterior = e.string_score;
      else if (t === "D-COLOC" && typeof e.pp4_min === "number") posterior = e.pp4_min;
      else if (t === "D-F2" && typeof e.loco_replication_fraction === "number")
        posterior = e.loco_replication_fraction;

      links.push({
        source: key,
        target: partner,
        layer,
        posterior,
        color: LAYER_COLOR_MAP[layer] ?? "#888",
        v4: v4Type
          ? {
              type: v4Type,
              emergence_stage: e.emergence_stage ?? e.emergence_stage_lr,
              loco_replication: e.loco_replication_fraction,
              pp4_min: e.pp4_min,
              druggable_pair: e.druggable_pair,
              conserved_mouse_a: e.conserved_mouse_a,
              conserved_mouse_b: e.conserved_mouse_b,
            }
          : undefined,
      });
    }
  }

  // Populate center node attributes from v2.attributes
  const center = nodeMap.get(key)!;
  center.biotype = v2.biotype ?? null;
  center.ensembl_id = v2.ensembl_id ?? "";
  center.is_deg = !!v2.attributes?.is_deg;
  center.is_conserved_core = !!v2.attributes?.is_conserved_core;
  center.dgidb_druggable = !!v2.attributes?.dgidb_druggable;
  center.dream_logfc = v2.attributes?.dream_logFC ?? null;
  center.dream_padj = v2.attributes?.dream_padj ?? null;
  center.sex_class = v2.attributes?.sex_class ?? null;
  center.degree = links.length;

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
