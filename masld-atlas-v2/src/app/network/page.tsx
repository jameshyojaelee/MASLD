"use client";

import { useEffect, useState, useCallback, useMemo, useRef, Suspense } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { ForceGraph } from "@/components/network/force-graph";
import { ControlsPanel, DEFAULT_THRESHOLDS, type Thresholds } from "@/components/network/controls-panel";
import { GeneDetailSidebar } from "@/components/network/gene-detail-sidebar";
import { FilterResults } from "@/components/network/filter-results";
import {
  executeFilter,
  filtersToParams,
  filtersFromParams,
  hydrateGeneSummaries,
  DEFAULT_FILTERS,
  filtersActive,
  type NetworkFilters,
  type FStageFilter,
  type GeneSummary,
} from "@/lib/network-filters";
import {
  loadSearchIndex,
  loadCommunities,
  loadGeneGraph,
  loadGeneGraphV2,
  loadGlobalLayout,
  loadLayerMetadata,
  loadCommunitiesF01,
  loadCommunitiesF34,
  DEFAULT_LAYERS,
} from "@/lib/network-data";
import type {
  CommunityData,
  CommunityInfo,
  EdgeLayer,
  EdgeType,
  GeneGraph,
  GeneGraphV2,
  LayerMetadata,
  NetworkEdge,
  NetworkNode,
  NetworkSearchEntry,
  NetworkViewMode,
  NodeColorBy,
} from "@/lib/network-types";

// Maximum queried genes. react-force-graph-2d begins to stutter beyond ~5 gene
// neighborhoods with >500 neighbors each; cap defensively.
const MAX_QUERIED_GENES = 5;

// ---------------------------------------------------------------------------
// Default layer set (all enabled)
// ---------------------------------------------------------------------------

const ALL_LAYERS: EdgeLayer[] = [
  "ppi", "coexpr", "regulon", "lr", "genetic",
  "pathway", "spatial", "cosmos", "cerna", "xspecies",
];

function parseLayersParam(param: string | null): Set<EdgeLayer> {
  if (!param) return new Set(ALL_LAYERS);
  const items = param.split(",").filter((s) => ALL_LAYERS.includes(s as EdgeLayer));
  return items.length > 0 ? new Set(items as EdgeLayer[]) : new Set(ALL_LAYERS);
}

function parseFdrParam(param: string | null): number {
  if (!param) return 0.1;
  const val = parseFloat(param);
  if (isNaN(val) || val < 0.01 || val > 0.5) return 0.1;
  return val;
}

/** Parse `genes=SYM1,SYM2,...` (preferred) or `gene=SYM` (single, legacy). */
function parseGenesParam(
  genes: string | null,
  gene: string | null
): string[] {
  const source = genes && genes.length > 0 ? genes : gene ?? "";
  if (!source) return [];
  return Array.from(
    new Set(
      source
        .split(",")
        .map((s) => s.trim().toUpperCase())
        .filter((s) => s.length > 0)
    )
  ).slice(0, MAX_QUERIED_GENES);
}

/** Parse the `th` URL param (compact, e.g. `0.70|0.70|0.50`). */
function parseThresholdsParam(param: string | null): Thresholds {
  if (!param) return DEFAULT_THRESHOLDS;
  const parts = param.split("|").map(parseFloat);
  if (parts.length !== 3 || parts.some((p) => isNaN(p))) return DEFAULT_THRESHOLDS;
  return { string: parts[0], loco: parts[1], pp4: parts[2] };
}

function serializeThresholds(t: Thresholds): string {
  return `${t.string.toFixed(2)}|${t.loco.toFixed(2)}|${t.pp4.toFixed(2)}`;
}

// ---------------------------------------------------------------------------
// View mode tabs
// ---------------------------------------------------------------------------

const VIEW_MODES: { key: NetworkViewMode; label: string; description: string }[] = [
  { key: "neighborhood", label: "Gene Neighborhood", description: "Explore edges around a single gene" },
  { key: "community", label: "Community Map", description: "Global layout colored by community" },
  { key: "compare", label: "Compare", description: "Overlay two gene neighborhoods" },
];

// ---------------------------------------------------------------------------
// Autocomplete search
// ---------------------------------------------------------------------------

function useGeneAutocomplete(searchIndex: NetworkSearchEntry[]) {
  const [query, setQuery] = useState("");
  const [suggestions, setSuggestions] = useState<NetworkSearchEntry[]>([]);

  useEffect(() => {
    if (query.length < 2 || searchIndex.length === 0) {
      setSuggestions([]);
      return;
    }

    const q = query.toUpperCase();
    const results = searchIndex
      .filter((g) => g.symbol.toUpperCase().includes(q))
      .sort((a, b) => {
        // Exact match first, then by degree
        const aExact = a.symbol.toUpperCase() === q ? 0 : 1;
        const bExact = b.symbol.toUpperCase() === q ? 0 : 1;
        if (aExact !== bExact) return aExact - bExact;
        return b.degree - a.degree;
      })
      .slice(0, 10);

    setSuggestions(results);
  }, [query, searchIndex]);

  return { query, setQuery, suggestions, setSuggestions };
}

// ---------------------------------------------------------------------------
// v2 filter state
// ---------------------------------------------------------------------------

const EDGE_TYPES_V2: EdgeType[] = ["S", "D-F2", "D-COLOC", "D-LR", "D-ceRNA", "D-XS"];

const F_STAGE_OPTIONS: { value: FStageFilter; label: string }[] = [
  { value: "all", label: "All" },
  { value: "F0-F1", label: "F0-F1" },
  { value: "F2", label: "F2" },
  { value: "F3-F4", label: "F3-F4" },
  { value: "F2_emerging", label: "F2 emerging" },
  { value: "F2_dissolving", label: "F2 dissolving" },
];

// ---------------------------------------------------------------------------
// Inner page component (needs Suspense boundary for useSearchParams)
// ---------------------------------------------------------------------------

function NetworkPageInner() {
  const searchParams = useSearchParams();
  const router = useRouter();

  // URL state
  const urlGenes = useMemo(
    () =>
      parseGenesParam(
        searchParams.get("genes"),
        searchParams.get("gene")
      ),
    [searchParams]
  );
  const urlGene = urlGenes[0] ?? null; // first queried gene = "primary" for v2 sidebar annotations
  const urlLayers = searchParams.get("layers");
  const urlTh = searchParams.get("th");
  const urlMode = searchParams.get("mode") as NetworkViewMode | null;
  const urlCompare = searchParams.get("compare") ?? null;
  const urlLabelTopK = searchParams.get("labelK");

  // Data state
  const [searchIndex, setSearchIndex] = useState<NetworkSearchEntry[]>([]);
  const [communities, setCommunities] = useState<CommunityData | null>(null);
  const [layerMeta, setLayerMeta] = useState<LayerMetadata[]>(DEFAULT_LAYERS);
  /**
   * Multi-gene graph state: one GeneGraph per queried symbol. Nodes/edges are
   * merged downstream via `useMemo` into a single graph. Preserves the v1-shape
   * `geneGraph` in the first slot for backcompat paths (e.g. the v2 sidebar).
   */
  const [geneGraphs, setGeneGraphs] = useState<Record<string, GeneGraph>>({});
  const [compareGraph, setCompareGraph] = useState<GeneGraph | null>(null);
  const [globalLayout, setGlobalLayout] = useState<NetworkNode[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // v2 filter + data state
  const [v2Filters, setV2Filters] = useState<NetworkFilters>(() =>
    filtersFromParams(new URLSearchParams(searchParams?.toString() ?? ""))
  );
  const [edgeTypesActive, setEdgeTypesActive] = useState<Set<EdgeType>>(
    new Set(EDGE_TYPES_V2)
  );
  const [filterResults, setFilterResults] = useState<GeneSummary[]>([]);
  const [filterUnresolved, setFilterUnresolved] = useState<string[]>([]);
  const [filterLoading, setFilterLoading] = useState(false);
  const filterSeqRef = useRef(0);
  const [v2Gene, setV2Gene] = useState<GeneGraphV2 | null>(null);
  const [communitiesF01, setCommunitiesF01] = useState<CommunityInfo[] | null>(null);
  const [communitiesF34, setCommunitiesF34] = useState<CommunityInfo[] | null>(null);
  const [showFilters, setShowFilters] = useState(true);

  // UI state
  const [viewMode, setViewMode] = useState<NetworkViewMode>(urlMode ?? "neighborhood");
  const [selectedLayers, setSelectedLayers] = useState<Set<EdgeLayer>>(
    parseLayersParam(urlLayers)
  );
  const [thresholds, setThresholds] = useState<Thresholds>(
    parseThresholdsParam(urlTh)
  );
  const [labelTopK, setLabelTopK] = useState<number>(() => {
    if (!urlLabelTopK) return 20;
    const v = parseInt(urlLabelTopK, 10);
    return isNaN(v) || v < 0 || v > 100 ? 20 : v;
  });
  const [colorBy, setColorBy] = useState<NodeColorBy>("community");
  const [selectedNode, setSelectedNode] = useState<NetworkNode | null>(null);
  const [, setHoveredNode] = useState<NetworkNode | null>(null);
  const [capError, setCapError] = useState<string | null>(null);

  // Search
  const { query, setQuery, suggestions, setSuggestions } = useGeneAutocomplete(searchIndex);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const searchRef = useRef<HTMLDivElement>(null);

  // Close suggestions on outside click
  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (searchRef.current && !searchRef.current.contains(e.target as Node)) {
        setShowSuggestions(false);
      }
    }
    document.addEventListener("mousedown", handleClick);
    return () => document.removeEventListener("mousedown", handleClick);
  }, []);

  // ---------------------------------------------------------------------------
  // Load static data on mount
  // ---------------------------------------------------------------------------

  useEffect(() => {
    loadSearchIndex().then(setSearchIndex).catch(() => {});
    loadCommunities().then(setCommunities).catch(() => {});
    loadLayerMetadata().then(setLayerMeta).catch(() => {});
    // v2 payloads (best-effort; network page gracefully degrades if missing)
    loadCommunitiesF01()
      .then((p) => setCommunitiesF01(p.communities))
      .catch(() => setCommunitiesF01(null));
    loadCommunitiesF34()
      .then((p) => setCommunitiesF34(p.communities))
      .catch(() => setCommunitiesF34(null));
  }, []);

  // Load v2 gene graph annotations when the URL gene changes
  useEffect(() => {
    if (!urlGene) {
      setV2Gene(null);
      return;
    }
    loadGeneGraphV2(urlGene)
      .then(setV2Gene)
      .catch(() => setV2Gene(null));
  }, [urlGene]);

  // ---------------------------------------------------------------------------
  // URL sync — update URL when state changes
  // ---------------------------------------------------------------------------

  const updateUrl = useCallback(
    (
      genes: string[],
      mode: NetworkViewMode,
      layers: Set<EdgeLayer>,
      th: Thresholds,
      compare: string | null,
      filters: NetworkFilters,
      topK: number
    ) => {
      const params = new URLSearchParams();
      if (genes.length === 1) {
        // Single gene: prefer `gene=` for backcompat with old share links.
        params.set("gene", genes[0]);
      } else if (genes.length > 1) {
        params.set("genes", genes.join(","));
      }
      if (mode !== "neighborhood") params.set("mode", mode);
      if (layers.size !== ALL_LAYERS.length) {
        params.set("layers", Array.from(layers).join(","));
      }
      const defaultTh =
        th.string === DEFAULT_THRESHOLDS.string &&
        th.loco === DEFAULT_THRESHOLDS.loco &&
        th.pp4 === DEFAULT_THRESHOLDS.pp4;
      if (!defaultTh) params.set("th", serializeThresholds(th));
      if (compare) params.set("compare", compare);
      if (topK !== 20) params.set("labelK", String(topK));
      for (const [k, v] of Object.entries(filtersToParams(filters))) {
        params.set(k, v);
      }
      const qs = params.toString();
      router.replace(qs ? `/network?${qs}` : "/network", { scroll: false });
    },
    [router]
  );

  // ---------------------------------------------------------------------------
  // Load gene graph when URL gene changes
  // ---------------------------------------------------------------------------

  const loadGenes = useCallback(
    async (symbols: string[]) => {
      if (symbols.length === 0) {
        setGeneGraphs({});
        return;
      }
      setLoading(true);
      setError(null);
      try {
        const results = await Promise.all(
          symbols.map(async (sym) => {
            try {
              return [sym, await loadGeneGraph(sym)] as const;
            } catch {
              return [sym, null] as const;
            }
          })
        );
        const next: Record<string, GeneGraph> = {};
        const missing: string[] = [];
        for (const [sym, g] of results) {
          if (g) next[sym] = g;
          else missing.push(sym);
        }
        setGeneGraphs(next);
        setSelectedNode(null);
        if (missing.length > 0) {
          setError(`Could not load network for: ${missing.join(", ")}`);
        }
      } finally {
        setLoading(false);
      }
    },
    []
  );

  useEffect(() => {
    loadGenes(urlGenes);
    // urlGenes is a memoized array — identity changes only when the actual
    // list changes, so this is safe.
  }, [urlGenes, loadGenes]);

  // Load compare gene if present
  useEffect(() => {
    if (urlCompare && viewMode === "compare") {
      loadGeneGraph(urlCompare)
        .then(setCompareGraph)
        .catch(() => setCompareGraph(null));
    } else {
      setCompareGraph(null);
    }
  }, [urlCompare, viewMode]);

  // Load global layout for community map mode
  useEffect(() => {
    if (viewMode === "community" && !globalLayout) {
      loadGlobalLayout().then(setGlobalLayout).catch(() => {});
    }
  }, [viewMode, globalLayout]);

  // ---------------------------------------------------------------------------
  // Handlers
  // ---------------------------------------------------------------------------

  /**
   * Replace the current query with a single gene (e.g. from autocomplete or
   * direct submit). For additive add, use `addGene`.
   */
  const selectGene = useCallback(
    (symbol: string) => {
      setQuery("");
      setShowSuggestions(false);
      setSuggestions([]);
      setCapError(null);
      updateUrl(
        [symbol.toUpperCase()],
        viewMode,
        selectedLayers,
        thresholds,
        urlCompare,
        v2Filters,
        labelTopK
      );
    },
    [viewMode, selectedLayers, thresholds, urlCompare, updateUrl, setQuery, setSuggestions, v2Filters, labelTopK]
  );

  /** Add a gene to the multi-gene query (chip input "Add gene" button). */
  const addGene = useCallback(
    (symbol: string) => {
      const sym = symbol.trim().toUpperCase();
      if (!sym) return;
      if (urlGenes.includes(sym)) {
        setQuery("");
        return;
      }
      if (urlGenes.length >= MAX_QUERIED_GENES) {
        setCapError(`Maximum ${MAX_QUERIED_GENES} genes at a time`);
        return;
      }
      setCapError(null);
      const next = [...urlGenes, sym];
      setQuery("");
      setShowSuggestions(false);
      setSuggestions([]);
      updateUrl(next, viewMode, selectedLayers, thresholds, urlCompare, v2Filters, labelTopK);
    },
    [urlGenes, viewMode, selectedLayers, thresholds, urlCompare, updateUrl, setQuery, setSuggestions, v2Filters, labelTopK]
  );

  const removeGene = useCallback(
    (symbol: string) => {
      const next = urlGenes.filter((g) => g !== symbol);
      setCapError(null);
      updateUrl(next, viewMode, selectedLayers, thresholds, urlCompare, v2Filters, labelTopK);
    },
    [urlGenes, viewMode, selectedLayers, thresholds, urlCompare, updateUrl, v2Filters, labelTopK]
  );

  const handleSearchSubmit = useCallback(
    (e: React.FormEvent) => {
      e.preventDefault();
      if (query.trim().length >= 2) {
        addGene(query.trim().toUpperCase());
      }
    },
    [query, addGene]
  );

  const handleToggleLayer = useCallback(
    (layer: EdgeLayer) => {
      setSelectedLayers((prev) => {
        const next = new Set(prev);
        if (next.has(layer)) {
          next.delete(layer);
        } else {
          next.add(layer);
        }
        updateUrl(urlGenes, viewMode, next, thresholds, urlCompare, v2Filters, labelTopK);
        return next;
      });
    },
    [urlGenes, viewMode, thresholds, urlCompare, updateUrl, v2Filters, labelTopK]
  );

  const handleThresholdsChange = useCallback(
    (value: Thresholds) => {
      setThresholds(value);
      updateUrl(urlGenes, viewMode, selectedLayers, value, urlCompare, v2Filters, labelTopK);
    },
    [urlGenes, viewMode, selectedLayers, urlCompare, updateUrl, v2Filters, labelTopK]
  );

  const handleLabelTopKChange = useCallback(
    (value: number) => {
      setLabelTopK(value);
      updateUrl(urlGenes, viewMode, selectedLayers, thresholds, urlCompare, v2Filters, value);
    },
    [urlGenes, viewMode, selectedLayers, thresholds, urlCompare, updateUrl, v2Filters]
  );

  const handleViewModeChange = useCallback(
    (mode: NetworkViewMode) => {
      setViewMode(mode);
      setSelectedNode(null);
      updateUrl(urlGenes, mode, selectedLayers, thresholds, urlCompare, v2Filters, labelTopK);
    },
    [urlGenes, selectedLayers, thresholds, urlCompare, updateUrl, v2Filters, labelTopK]
  );

  const handleNodeClick = useCallback((node: NetworkNode) => {
    setSelectedNode(node);
  }, []);

  const handleNodeHover = useCallback((node: NetworkNode | null) => {
    setHoveredNode(node);
  }, []);

  const handleCloseSidebar = useCallback(() => {
    setSelectedNode(null);
  }, []);

  const updateFilters = useCallback(
    (updater: (prev: NetworkFilters) => NetworkFilters) => {
      setV2Filters((prev) => {
        const next = updater(prev);
        updateUrl(urlGenes, viewMode, selectedLayers, thresholds, urlCompare, next, labelTopK);
        return next;
      });
    },
    [urlGenes, viewMode, selectedLayers, thresholds, urlCompare, updateUrl, labelTopK]
  );

  const clearFilters = useCallback(() => {
    updateFilters(() => DEFAULT_FILTERS);
  }, [updateFilters]);

  // Execute filter whenever v2Filters changes (debounced via seq counter).
  // Retries automatically if filter indexes not yet deployed.
  useEffect(() => {
    const seq = ++filterSeqRef.current;
    let cancelled = false;
    let retryTimer: ReturnType<typeof setTimeout> | null = null;

    async function run() {
      if (!filtersActive(v2Filters)) {
        setFilterResults([]);
        setFilterUnresolved([]);
        setFilterLoading(false);
        return;
      }
      setFilterLoading(true);
      try {
        const { genes, unresolved } = await executeFilter(v2Filters);
        if (cancelled || seq !== filterSeqRef.current) return;
        const summaries = await hydrateGeneSummaries(genes);
        if (cancelled || seq !== filterSeqRef.current) return;
        summaries.sort((a, b) => (b.degree ?? 0) - (a.degree ?? 0));
        setFilterResults(summaries);
        setFilterUnresolved(unresolved);
        setFilterLoading(false);
        if (unresolved.length > 0) {
          // Retry in 5s in case indexes are still being deployed
          retryTimer = setTimeout(() => {
            if (!cancelled && seq === filterSeqRef.current) run();
          }, 5000);
        }
      } catch {
        if (!cancelled && seq === filterSeqRef.current) {
          setFilterLoading(false);
        }
      }
    }
    run();
    return () => {
      cancelled = true;
      if (retryTimer) clearTimeout(retryTimer);
    };
  }, [v2Filters]);

  // ---------------------------------------------------------------------------
  // Merge graph data based on view mode
  // ---------------------------------------------------------------------------

  const queriedIdSet = useMemo(() => new Set(urlGenes), [urlGenes]);

  /**
   * Merge graph data across all queried gene neighborhoods. For each non-queried
   * neighbor, attach `sharedCount` = number of queried genes that have it as a
   * neighbor. Queried nodes get `queried = true`. Edges are deduped by
   * `source|target|layer`.
   */
  const { nodes, links, summary } = useMemo(() => {
    if (viewMode === "community" && globalLayout) {
      return {
        nodes: globalLayout,
        links: [] as NetworkEdge[],
        summary: null as null | {
          queried: number;
          shared: number;
          bridges: number;
          directEdges: number;
        },
      };
    }

    const graphs = urlGenes.map((g) => geneGraphs[g]).filter(Boolean) as GeneGraph[];

    if (viewMode === "compare" && graphs[0] && compareGraph) {
      const nodeMap = new Map<string, NetworkNode>();
      for (const n of graphs[0].nodes) nodeMap.set(n.id, n);
      for (const n of compareGraph.nodes) {
        if (!nodeMap.has(n.id)) nodeMap.set(n.id, n);
      }
      const mergedLinks = [...graphs[0].links, ...compareGraph.links];
      const linkSet = new Set<string>();
      const dedupLinks = mergedLinks.filter((l) => {
        const key = [l.source, l.target, l.layer].sort().join("|");
        if (linkSet.has(key)) return false;
        linkSet.add(key);
        return true;
      });
      return { nodes: Array.from(nodeMap.values()), links: dedupLinks, summary: null };
    }

    if (graphs.length === 0) {
      return { nodes: [] as NetworkNode[], links: [] as NetworkEdge[], summary: null };
    }

    // Count neighborhood membership across queried genes
    const neighborOwners = new Map<string, Set<string>>();
    for (const g of graphs) {
      const center = g.center;
      for (const n of g.nodes) {
        if (n.id === center) continue;
        if (!neighborOwners.has(n.id)) neighborOwners.set(n.id, new Set());
        neighborOwners.get(n.id)!.add(center);
      }
    }

    const nodeMap = new Map<string, NetworkNode & { sharedCount?: number; queried?: boolean }>();
    for (const g of graphs) {
      for (const n of g.nodes) {
        if (nodeMap.has(n.id)) continue;
        const queried = queriedIdSet.has(n.id);
        const sharedCount = queried ? urlGenes.length : neighborOwners.get(n.id)?.size ?? 0;
        nodeMap.set(n.id, { ...n, queried, sharedCount });
      }
    }

    // Dedup links across graphs
    const linkSet = new Set<string>();
    const mergedLinks: NetworkEdge[] = [];
    for (const g of graphs) {
      for (const l of g.links) {
        const srcId = typeof l.source === "object" ? (l.source as NetworkNode).id : l.source;
        const tgtId = typeof l.target === "object" ? (l.target as NetworkNode).id : l.target;
        const key = [srcId, tgtId, l.layer].sort().join("|");
        if (linkSet.has(key)) continue;
        linkSet.add(key);
        mergedLinks.push(l);
      }
    }

    const sharedCount = Array.from(neighborOwners.values()).filter(
      (s) => s.size >= 2
    ).length;
    const directEdges = mergedLinks.filter((l) => {
      const srcId = typeof l.source === "object" ? (l.source as NetworkNode).id : l.source;
      const tgtId = typeof l.target === "object" ? (l.target as NetworkNode).id : l.target;
      return queriedIdSet.has(srcId) && queriedIdSet.has(tgtId);
    }).length;

    return {
      nodes: Array.from(nodeMap.values()),
      links: mergedLinks,
      summary: {
        queried: urlGenes.length,
        shared: sharedCount,
        bridges: sharedCount, // bridge = shared (alias: neighbors of >=2 queried)
        directEdges,
      },
    };
  }, [viewMode, urlGenes, geneGraphs, compareGraph, globalLayout, queriedIdSet]);

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  return (
    <div className="flex h-screen flex-col overflow-hidden">
      {/* Top bar */}
      <div className="shrink-0 space-y-3 border-b border-border px-6 py-4">
        {/* Title + search + view tabs */}
        <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
          <div className="flex items-center gap-4">
            <h1 className="text-2xl font-bold tracking-tight">Network Explorer</h1>

            {/* Gene chip input: add-genes form */}
            <div ref={searchRef} className="relative w-80">
              <form onSubmit={handleSearchSubmit} className="flex gap-1">
                <Input
                  placeholder={
                    urlGenes.length === 0
                      ? "Search gene (e.g. THRB)..."
                      : `Add gene (${urlGenes.length}/${MAX_QUERIED_GENES})...`
                  }
                  value={query}
                  onChange={(e) => {
                    setQuery(e.target.value);
                    setShowSuggestions(true);
                  }}
                  onFocus={() => setShowSuggestions(true)}
                  className="flex-1"
                />
                <Button
                  type="submit"
                  size="sm"
                  variant="secondary"
                  disabled={query.trim().length < 2 || urlGenes.length >= MAX_QUERIED_GENES}
                >
                  Add gene
                </Button>
              </form>

              {/* Autocomplete dropdown */}
              {showSuggestions && suggestions.length > 0 && (
                <div className="absolute left-0 top-full z-50 mt-1 w-full overflow-hidden rounded-lg border border-border bg-popover shadow-lg">
                  <ul className="max-h-64 overflow-y-auto py-1">
                    {suggestions.map((s) => (
                      <li key={s.symbol}>
                        <button
                          type="button"
                          onClick={() => addGene(s.symbol)}
                          className="flex w-full items-center justify-between px-3 py-1.5 text-sm transition-colors hover:bg-muted"
                        >
                          <span className="font-mono font-medium">{s.symbol}</span>
                          <span className="text-xs text-muted-foreground">
                            deg: {s.degree}
                            {s.community_label && ` | ${s.community_label}`}
                          </span>
                        </button>
                      </li>
                    ))}
                  </ul>
                </div>
              )}
            </div>
          </div>

          {/* View mode tabs */}
          <div className="flex gap-1">
            {VIEW_MODES.map((mode) => (
              <Button
                key={mode.key}
                variant={viewMode === mode.key ? "default" : "outline"}
                size="sm"
                onClick={() => handleViewModeChange(mode.key)}
                title={mode.description}
              >
                {mode.label}
              </Button>
            ))}
          </div>
        </div>

        {/* Controls panel */}
        <ControlsPanel
          layers={layerMeta}
          selectedLayers={selectedLayers}
          onToggleLayer={handleToggleLayer}
          thresholds={thresholds}
          onThresholdsChange={handleThresholdsChange}
          colorBy={colorBy}
          onColorByChange={setColorBy}
          labelTopK={labelTopK}
          onLabelTopKChange={handleLabelTopKChange}
        />

        {/* Queried-genes chip strip */}
        {urlGenes.length > 0 && (
          <div className="flex flex-wrap items-center gap-2 text-sm">
            <span className="text-muted-foreground">Querying:</span>
            {urlGenes.map((g) => (
              <Badge key={g} variant="secondary" className="gap-1 font-mono">
                {g}
                <button
                  type="button"
                  onClick={() => removeGene(g)}
                  aria-label={`Remove ${g}`}
                  className="ml-0.5 text-muted-foreground hover:text-foreground"
                >
                  &#x2715;
                </button>
              </Badge>
            ))}
            {viewMode === "compare" && urlCompare && (
              <>
                <span className="text-muted-foreground">vs</span>
                <Badge variant="secondary" className="font-mono">
                  {urlCompare}
                </Badge>
              </>
            )}
            {loading && (
              <span className="text-xs text-muted-foreground animate-pulse">Loading...</span>
            )}
            {error && <span className="text-xs text-destructive">{error}</span>}
          </div>
        )}
        {capError && (
          <div className="rounded-md border border-destructive bg-destructive/10 px-3 py-1.5 text-xs text-destructive">
            {capError}
          </div>
        )}
      </div>

      {/* Main content: filter sidebar + graph + detail sidebar */}
      <div className="flex min-h-0 flex-1">
        {/* v2 filter sidebar */}
        {showFilters && (
          <aside className="w-64 shrink-0 overflow-y-auto border-r border-border bg-muted/30 p-4 text-sm">
            <div className="mb-3 flex items-center justify-between">
              <h2 className="font-semibold">Filters</h2>
              <button
                type="button"
                onClick={() => setShowFilters(false)}
                className="text-xs text-muted-foreground hover:text-foreground"
              >
                hide
              </button>
            </div>

            <fieldset className="mb-4">
              <legend className="mb-1 text-xs font-medium uppercase text-muted-foreground">F-stage</legend>
              <div className="flex flex-wrap gap-1">
                {F_STAGE_OPTIONS.map((opt) => (
                  <Button
                    key={opt.value}
                    size="sm"
                    variant={v2Filters.fstage === opt.value ? "default" : "outline"}
                    onClick={() => updateFilters((p) => ({ ...p, fstage: opt.value }))}
                  >
                    {opt.label}
                  </Button>
                ))}
              </div>
            </fieldset>

            <fieldset className="mb-4">
              <legend className="mb-1 text-xs font-medium uppercase text-muted-foreground">Edge type (display)</legend>
              <div className="flex flex-col gap-1">
                {EDGE_TYPES_V2.map((t) => {
                  const checked = edgeTypesActive.has(t);
                  return (
                    <label key={t} className="flex cursor-pointer items-center gap-2">
                      <input
                        type="checkbox"
                        checked={checked}
                        onChange={() =>
                          setEdgeTypesActive((prev) => {
                            const next = new Set(prev);
                            if (next.has(t)) next.delete(t);
                            else next.add(t);
                            return next;
                          })
                        }
                      />
                      <span className="font-mono">{t}</span>
                    </label>
                  );
                })}
              </div>
            </fieldset>

            <fieldset className="mb-4 flex flex-col gap-2">
              <legend className="mb-1 text-xs font-medium uppercase text-muted-foreground">Flags</legend>
              <label className="flex cursor-pointer items-center gap-2">
                <input
                  type="checkbox"
                  checked={v2Filters.hasColoc}
                  onChange={(e) => updateFilters((p) => ({ ...p, hasColoc: e.target.checked }))}
                />
                <span>Has COLOC (PP4 &ge; 0.5)</span>
              </label>
              <label className="flex cursor-pointer items-center gap-2">
                <input
                  type="checkbox"
                  checked={v2Filters.druggable}
                  onChange={(e) => updateFilters((p) => ({ ...p, druggable: e.target.checked }))}
                />
                <span>Druggable gene</span>
              </label>
              <label className="flex cursor-pointer items-center gap-2">
                <input
                  type="checkbox"
                  checked={v2Filters.druggablePair}
                  onChange={(e) => updateFilters((p) => ({ ...p, druggablePair: e.target.checked }))}
                />
                <span>Druggable pair (both endpoints)</span>
              </label>
            </fieldset>

            <fieldset className="mb-4">
              <legend className="mb-1 text-xs font-medium uppercase text-muted-foreground">
                Cohort support &ge; {v2Filters.cohortSupport}
              </legend>
              <input
                type="range"
                min={0}
                max={10}
                step={1}
                value={v2Filters.cohortSupport}
                onChange={(e) =>
                  updateFilters((p) => ({ ...p, cohortSupport: parseInt(e.target.value, 10) }))
                }
                className="w-full"
              />
            </fieldset>

            {(communitiesF01 || communitiesF34) && (
              <fieldset className="mb-4 flex flex-col gap-2">
                <legend className="mb-1 text-xs font-medium uppercase text-muted-foreground">
                  Community
                </legend>
                {communitiesF01 && (
                  <label className="flex items-center gap-2 text-xs">
                    <span className="w-14 text-muted-foreground">F0-F1</span>
                    <select
                      value={v2Filters.communityF01 ?? ""}
                      onChange={(e) =>
                        updateFilters((p) => ({
                          ...p,
                          communityF01: e.target.value === "" ? null : parseInt(e.target.value, 10),
                        }))
                      }
                      className="h-7 flex-1 rounded-md border border-input bg-transparent px-1 text-xs"
                    >
                      <option value="">Any</option>
                      {communitiesF01.slice(0, 10).map((c) => (
                        <option key={c.id} value={c.id}>
                          C{c.id} ({c.n_genes})
                          {c.top_hallmark ? " - " + c.top_hallmark.replace(/^HALLMARK_/, "").slice(0, 20) : ""}
                        </option>
                      ))}
                    </select>
                  </label>
                )}
                {communitiesF34 && (
                  <label className="flex items-center gap-2 text-xs">
                    <span className="w-14 text-muted-foreground">F3-F4</span>
                    <select
                      value={v2Filters.communityF34 ?? ""}
                      onChange={(e) =>
                        updateFilters((p) => ({
                          ...p,
                          communityF34: e.target.value === "" ? null : parseInt(e.target.value, 10),
                        }))
                      }
                      className="h-7 flex-1 rounded-md border border-input bg-transparent px-1 text-xs"
                    >
                      <option value="">Any</option>
                      {communitiesF34.slice(0, 10).map((c) => (
                        <option key={c.id} value={c.id}>
                          C{c.id} ({c.n_genes})
                          {c.top_hallmark ? " - " + c.top_hallmark.replace(/^HALLMARK_/, "").slice(0, 20) : ""}
                        </option>
                      ))}
                    </select>
                  </label>
                )}
              </fieldset>
            )}

            <FilterResults
              loading={filterLoading}
              results={filterResults}
              unresolved={filterUnresolved}
              selectedGene={urlGene}
              onSelect={selectGene}
              onClear={clearFilters}
            />

            {/* Community membership summary */}
            {v2Gene?.community && (
              <section className="mt-4 border-t border-border pt-3">
                <h3 className="mb-2 text-xs font-medium uppercase text-muted-foreground">
                  Community
                </h3>
                {v2Gene.community.f01 && (
                  <div className="mb-2">
                    <Badge variant="outline" className="mr-1">F0-F1</Badge>
                    <span className="font-mono text-xs">
                      C{v2Gene.community.f01.macro_id ?? v2Gene.community.f01.community_id}
                    </span>
                    {v2Gene.community.f01.top_hallmark && (
                      <div className="mt-0.5 text-xs text-muted-foreground">
                        {v2Gene.community.f01.top_hallmark.replace(/^HALLMARK_/, "")}
                      </div>
                    )}
                  </div>
                )}
                {v2Gene.community.f34 && (
                  <div>
                    <Badge variant="outline" className="mr-1">F3-F4</Badge>
                    <span className="font-mono text-xs">
                      C{v2Gene.community.f34.macro_id ?? v2Gene.community.f34.community_id}
                    </span>
                    {v2Gene.community.f34.top_hallmark && (
                      <div className="mt-0.5 text-xs text-muted-foreground">
                        {v2Gene.community.f34.top_hallmark.replace(/^HALLMARK_/, "")}
                      </div>
                    )}
                  </div>
                )}
              </section>
            )}

            {/* Community lists (counts) */}
            {(communitiesF01 || communitiesF34) && (
              <section className="mt-4 border-t border-border pt-3 text-xs text-muted-foreground">
                {communitiesF01 && <div>F0-F1 communities: {communitiesF01.length}</div>}
                {communitiesF34 && <div>F3-F4 communities: {communitiesF34.length}</div>}
              </section>
            )}
          </aside>
        )}
        {!showFilters && (
          <button
            type="button"
            onClick={() => setShowFilters(true)}
            className="absolute left-2 top-32 z-20 rounded-md border border-border bg-background px-2 py-1 text-xs shadow hover:bg-muted"
          >
            Show filters
          </button>
        )}

        {/* Graph canvas */}
        <div className="min-w-0 flex-1 flex flex-col">
          {/* Multi-gene summary bar */}
          {summary && urlGenes.length > 0 && (
            <div className="shrink-0 border-b border-border bg-muted/30 px-4 py-1.5 text-xs text-muted-foreground">
              <span className="font-medium text-foreground">{summary.queried}</span> queried gene
              {summary.queried === 1 ? "" : "s"} &middot;{" "}
              <span className="font-medium text-foreground">{summary.shared}</span> shared
              neighbors (&ge;2 queries) &middot;{" "}
              <span className="font-medium text-foreground">{summary.bridges}</span> bridge genes
              &middot;{" "}
              <span className="font-medium text-foreground">{summary.directEdges}</span> direct
              edges
            </div>
          )}
          <div className="min-h-0 flex-1">
            <ForceGraph
              nodes={nodes}
              links={links}
              onNodeClick={handleNodeClick}
              onNodeHover={handleNodeHover}
              selectedLayers={selectedLayers}
              thresholds={thresholds}
              colorBy={colorBy}
              viewMode={viewMode}
              highlightNodeId={urlGene ?? undefined}
              queriedIds={queriedIdSet}
              labelTopK={labelTopK}
              fstage={v2Filters.fstage}
            />
          </div>
        </div>

        {/* Gene detail sidebar (slides in when a node is selected) */}
        {selectedNode && (
          <GeneDetailSidebar
            node={selectedNode}
            neighbors={nodes}
            edges={links}
            communities={communities}
            onClose={handleCloseSidebar}
          />
        )}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page wrapper with Suspense for useSearchParams
// ---------------------------------------------------------------------------

/** Layout-shaped skeleton (top bar + filter rail + graph viewport). */
function NetworkSkeleton() {
  return (
    <div className="flex h-screen flex-col overflow-hidden">
      <div className="shrink-0 space-y-3 border-b border-border px-6 py-4">
        <div className="flex items-center gap-4">
          <div className="h-7 w-44 animate-pulse rounded bg-muted" />
          <div className="h-9 w-80 animate-pulse rounded bg-muted" />
        </div>
        <div className="h-8 w-full max-w-2xl animate-pulse rounded bg-muted/70" />
      </div>
      <div className="flex min-h-0 flex-1">
        <div className="hidden w-64 shrink-0 space-y-3 border-r border-border bg-muted/20 p-4 sm:block">
          {Array.from({ length: 5 }).map((_, i) => (
            <div key={i} className="h-8 w-full animate-pulse rounded bg-muted" />
          ))}
        </div>
        <div className="relative min-w-0 flex-1 overflow-hidden bg-card">
          <svg
            className="h-full w-full animate-pulse text-muted-foreground/40"
            viewBox="0 0 100 100"
            preserveAspectRatio="xMidYMid slice"
            aria-hidden
          >
            {[
              [35, 30],
              [66, 28],
              [30, 62],
              [70, 66],
              [50, 72],
              [20, 46],
              [82, 48],
            ].map(([x, y], i) => (
              <line key={i} x1={50} y1={45} x2={x} y2={y} stroke="currentColor" strokeWidth={0.4} />
            ))}
            {[
              [50, 45, 5],
              [35, 30, 3],
              [66, 28, 3.5],
              [30, 62, 3],
              [70, 66, 4],
              [50, 72, 2.5],
              [20, 46, 2.5],
              [82, 48, 2.5],
            ].map(([x, y, r], i) => (
              <circle key={i} cx={x} cy={y} r={r} fill="currentColor" />
            ))}
          </svg>
          <span className="absolute bottom-3 left-3 text-sm text-muted-foreground">
            Loading network explorer…
          </span>
        </div>
      </div>
    </div>
  );
}

export default function NetworkPage() {
  return (
    <Suspense fallback={<NetworkSkeleton />}>
      <NetworkPageInner />
    </Suspense>
  );
}
