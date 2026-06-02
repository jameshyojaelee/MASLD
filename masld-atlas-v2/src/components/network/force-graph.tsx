"use client";

/**
 * Force-directed graph wrapper around react-force-graph-2d.
 *
 * NOTE: `react-force-graph` must be added to package.json dependencies:
 *   "react-force-graph": "^1.46.0"
 * Then run: npm install
 *
 * The package is dynamically imported below to avoid SSR issues (it requires
 * a browser canvas context). The component renders a loading skeleton until
 * the library is available.
 */

import {
  useEffect,
  useRef,
  useState,
  useCallback,
  useMemo,
} from "react";
import dynamic from "next/dynamic";
import { getLayerColor } from "@/lib/network-data";
import type { Thresholds } from "@/components/network/controls-panel";
import type {
  EdgeLayer,
  NetworkNode,
  NetworkEdge,
  NetworkViewMode,
  NodeColorBy,
} from "@/lib/network-types";

// Dynamic import to bypass SSR — react-force-graph requires <canvas>
const ForceGraph2D = dynamic(() => import("react-force-graph-2d"), {
  ssr: false,
  loading: () => <GraphSkeleton />,
});

// ---------------------------------------------------------------------------
// Props
// ---------------------------------------------------------------------------

interface ForceGraphProps {
  nodes: NetworkNode[];
  links: NetworkEdge[];
  onNodeClick: (node: NetworkNode) => void;
  onNodeHover: (node: NetworkNode | null) => void;
  selectedLayers: Set<EdgeLayer>;
  thresholds: Thresholds;
  colorBy: NodeColorBy;
  viewMode: NetworkViewMode;
  highlightNodeId?: string | null;
  /**
   * Set of gene symbols in the current multi-gene query. Queried nodes render
   * with a bold outline and are always labeled. Nodes whose `sharedCount >= 2`
   * (i.e. shared neighbors of >=2 queried genes) render at full opacity.
   */
  queriedIds?: Set<string>;
  /**
   * Top-K node-label budget (sorted by degree). 0 disables top-K labeling —
   * only queried/hovered nodes will show labels. Default 20.
   */
  labelTopK?: number;
  /**
   * Active F-stage filter from the sidebar (e.g. "F2_emerging"). When set to
   * anything other than "all", D-F2 edges whose emergence_stage does not match
   * are hidden (and orphan nodes are pruned so the graph doesn't disperse).
   */
  fstage?: string;
  width?: number;
  height?: number;
}

/**
 * Map the F-stage sidebar filter values onto the emergence_stage field
 * present on v4 D-F2 edges. Returns the set of emergence_stage strings that
 * should survive when the given filter is active, or `null` to disable the
 * filter entirely.
 */
function fstageEdgeFilter(fstage: string | undefined): Set<string> | null {
  if (!fstage || fstage === "all") return null;
  switch (fstage) {
    case "F2_emerging":
      return new Set(["F2_emerging"]);
    case "F2_dissolving":
      return new Set(["F2_dissolving"]);
    case "F0-F1":
      // F0-F1 neighborhoods: keep F0-specific + dissolving (stronger in F0-F1)
      return new Set(["F0_specific", "F2_dissolving"]);
    case "F2":
      return new Set(["F2_emerging", "F2_dissolving", "transient_F2"]);
    case "F3-F4":
      return new Set(["F34_specific", "F2_emerging", "progressive_up", "progressive_down"]);
    default:
      return null;
  }
}

/**
 * Per-edge-type visibility dispatch. v4 edges carry a `type` field that names
 * which threshold applies. Non-v4 edges fall through at full opacity.
 *
 * Returns `show: true` always so users retain context — edges below the
 * threshold are `dim: true` and get rendered at alpha 0.15.
 */
function edgeVisible(
  link: NetworkEdge,
  th: Thresholds
): { show: boolean; dim: boolean } {
  const t = link.v4?.type;
  if (t === "S") {
    const ok = link.posterior >= th.string;
    return { show: true, dim: !ok };
  }
  if (t === "D-F2") {
    const loco = link.v4?.loco_replication ?? 0;
    return { show: true, dim: loco < th.loco };
  }
  if (t === "D-COLOC") {
    const pp4 = link.v4?.pp4_min ?? 0;
    return { show: true, dim: pp4 < th.pp4 };
  }
  return { show: true, dim: false };
}

// ---------------------------------------------------------------------------
// Color helpers
// ---------------------------------------------------------------------------

const COMMUNITY_PALETTE = [
  "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6",
  "#1abc9c", "#e67e22", "#34495e", "#e91e63", "#607d8b",
  "#ff6b6b", "#4ecdc4", "#45b7d1", "#96ceb4", "#ffeaa7",
];

const SEX_CLASS_COLORS: Record<string, string> = {
  Female_biased: "#e74c3c",
  Male_biased: "#3498db",
  Divergent: "#9b59b6",
  Concordant: "#95a5a6",
};

const PROGRESSION_COLORS: Record<string, string> = {
  onset: "#f39c12",
  progression: "#e74c3c",
  late: "#9b59b6",
};

function getNodeColor(node: NetworkNode, colorBy: NodeColorBy): string {
  switch (colorBy) {
    case "community":
      return COMMUNITY_PALETTE[node.community_macro % COMMUNITY_PALETTE.length];
    case "layers":
      // Gradient from gray (0) to primary blue (7)
      return `hsl(220, ${Math.round((node.layers_active / 7) * 80)}%, ${60 - Math.round((node.layers_active / 7) * 20)}%)`;
    case "sex_class":
      return SEX_CLASS_COLORS[node.sex_class ?? ""] ?? "#95a5a6";
    case "progression":
      return PROGRESSION_COLORS[node.progression_class ?? ""] ?? "#95a5a6";
    case "logfc": {
      if (node.dream_logfc == null) return "#95a5a6";
      if (node.dream_logfc > 0) {
        const intensity = Math.min(1, Math.abs(node.dream_logfc) / 2);
        return `hsl(0, ${Math.round(intensity * 80)}%, ${60 - Math.round(intensity * 15)}%)`;
      }
      const intensity = Math.min(1, Math.abs(node.dream_logfc) / 2);
      return `hsl(220, ${Math.round(intensity * 80)}%, ${60 - Math.round(intensity * 15)}%)`;
    }
    case "druggability":
      return node.dgidb_druggable ? "#2ecc71" : "#95a5a6";
    default:
      return "#95a5a6";
  }
}

function getNodeSize(node: NetworkNode): number {
  // Base size 2, gentle log-scale on degree — keeps large hubs readable
  // without overwhelming the canvas.
  return 2 + Math.log2(Math.max(1, node.degree)) * 0.9;
}

// ---------------------------------------------------------------------------
// Skeleton placeholder
// ---------------------------------------------------------------------------

function GraphSkeleton() {
  return (
    <div className="flex h-full w-full items-center justify-center bg-card">
      <div className="flex flex-col items-center gap-3 text-muted-foreground">
        <svg className="size-8 animate-spin" fill="none" viewBox="0 0 24 24">
          <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
          <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z" />
        </svg>
        <span className="text-sm">Loading network graph...</span>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

export function ForceGraph({
  nodes,
  links,
  onNodeClick,
  onNodeHover,
  selectedLayers,
  thresholds,
  colorBy,
  viewMode,
  highlightNodeId,
  queriedIds,
  labelTopK = 20,
  fstage,
  width,
  height,
}: ForceGraphProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const fgRef = useRef<any>(null);
  const [hoveredNode, setHoveredNode] = useState<string | null>(null);
  const [dimensions, setDimensions] = useState({ w: width ?? 800, h: height ?? 600 });

  // Measure container
  useEffect(() => {
    if (width && height) {
      setDimensions({ w: width, h: height });
      return;
    }

    const el = containerRef.current;
    if (!el) return;

    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width: w, height: h } = entry.contentRect;
        if (w > 0 && h > 0) {
          setDimensions({ w: Math.floor(w), h: Math.floor(h) });
        }
      }
    });

    observer.observe(el);
    // Initial measurement
    const rect = el.getBoundingClientRect();
    if (rect.width > 0 && rect.height > 0) {
      setDimensions({ w: Math.floor(rect.width), h: Math.floor(rect.height) });
    }

    return () => observer.disconnect();
  }, [width, height]);

  // Filter links by selected edge layer + F-stage filter. Threshold dispatch
  // is handled in paintLink via edgeVisible() — below-threshold edges dim
  // rather than hide. F-stage applies only to D-F2 edges (other edge types
  // pass through unaffected regardless of F-stage selection).
  const fstageSet = useMemo(() => fstageEdgeFilter(fstage), [fstage]);

  // In multi-gene mode we keep only queried genes + shared neighbors (seen in
  // >=2 queried-gene neighborhoods). This anchors the graph into a single
  // connected core instead of letting pairwise-specific neighbors separate
  // into disconnected clusters that fly apart under repulsion.
  const multiGeneMode = !!queriedIds && queriedIds.size >= 2;
  const keepSharedOnly = multiGeneMode;

  const sharedKeepIds = useMemo(() => {
    if (!keepSharedOnly) return null;
    const keep = new Set<string>();
    if (queriedIds) for (const q of queriedIds) keep.add(q);
    for (const n of nodes) {
      const sc = (n as unknown as { sharedCount?: number }).sharedCount ?? 0;
      if (sc >= 2) keep.add(n.id);
    }
    return keep;
  }, [keepSharedOnly, queriedIds, nodes]);

  const filteredLinks = useMemo(() => {
    return links.filter((link) => {
      if (!selectedLayers.has(link.layer)) return false;
      if (fstageSet && link.v4?.type === "D-F2") {
        const es = link.v4.emergence_stage;
        if (!es || !fstageSet.has(es)) return false;
      }
      if (sharedKeepIds) {
        const src = typeof link.source === "object" ? (link.source as NetworkNode).id : link.source;
        const tgt = typeof link.target === "object" ? (link.target as NetworkNode).id : link.target;
        if (!sharedKeepIds.has(src) || !sharedKeepIds.has(tgt)) return false;
      }
      return true;
    });
  }, [links, selectedLayers, fstageSet, sharedKeepIds]);

  // Set of visible node IDs = endpoints of visible edges. The center/queried
  // genes are only kept if they participate in at least one surviving edge —
  // dropping them otherwise prevents the "orphan drift" where a filter removes
  // all of a gene's edges but the node stays, floating untethered. In
  // multi-gene mode queried genes are always kept (even if no shared-neighbor
  // edges survive) so the user can see which queries are isolated.
  const visibleNodeIds = useMemo(() => {
    const ids = new Set<string>();
    for (const link of filteredLinks) {
      ids.add(typeof link.source === "object" ? (link.source as NetworkNode).id : link.source);
      ids.add(typeof link.target === "object" ? (link.target as NetworkNode).id : link.target);
    }
    if (multiGeneMode && queriedIds) {
      for (const q of queriedIds) ids.add(q);
    }
    return ids;
  }, [filteredLinks, multiGeneMode, queriedIds]);

  const filteredNodes = useMemo(() => {
    if (viewMode === "community") return nodes;
    return nodes.filter((n) => visibleNodeIds.has(n.id));
  }, [nodes, visibleNodeIds, viewMode]);

  // Build graph data object for ForceGraph2D
  const graphData = useMemo(
    () => ({ nodes: filteredNodes, links: filteredLinks }),
    [filteredNodes, filteredLinks]
  );

  /**
   * Node-label policy: always label (a) queried genes, (b) highlight/center,
   * (c) top-K nodes by degree. Hover still reveals labels for any node via
   * the per-frame `showLabel` branch in paintNode.
   */
  const alwaysLabeledIds = useMemo(() => {
    const ids = new Set<string>();
    if (highlightNodeId) ids.add(highlightNodeId);
    if (queriedIds) for (const q of queriedIds) ids.add(q);
    if (labelTopK > 0) {
      const sorted = [...filteredNodes]
        .sort((a, b) => (b.degree ?? 0) - (a.degree ?? 0))
        .slice(0, labelTopK);
      for (const n of sorted) ids.add(n.id);
    }
    return ids;
  }, [filteredNodes, labelTopK, highlightNodeId, queriedIds]);

  // Precompute neighbor-sharedCount lookups. `sharedCount` is attached to the
  // merged node by the parent page during multi-gene merge; we read it here
  // with a safe fallback so this component stays independent of that merge.
  const sharedCountOf = useCallback(
    (n: NetworkNode): number =>
      (n as unknown as { sharedCount?: number }).sharedCount ?? 0,
    []
  );

  // Handle hover
  const handleNodeHover = useCallback(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (node: any) => {
      setHoveredNode(node?.id ?? null);
      onNodeHover(node as NetworkNode | null);
    },
    [onNodeHover]
  );

  // Handle click
  const handleNodeClick = useCallback(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (node: any) => {
      onNodeClick(node as NetworkNode);
    },
    [onNodeClick]
  );

  // Node rendering
  const paintNode = useCallback(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (node: any, ctx: CanvasRenderingContext2D, globalScale: number) => {
      const n = node as NetworkNode;
      const r = getNodeSize(n);
      const isHighlighted = n.id === highlightNodeId;
      const isHovered = n.id === hoveredNode;
      const isQueried = !!queriedIds && queriedIds.has(n.id);
      const shared = sharedCountOf(n);

      // Opacity: queried + shared (>=2) nodes at full opacity; queried singles
      // also full; otherwise fade to 0.4 when the user is running a multi-gene
      // query (so shared neighbors pop). If no queriedIds, default behavior.
      const hasQuery = !!queriedIds && queriedIds.size > 0;
      const nodeAlpha = !hasQuery || isQueried || shared >= 2 ? 1.0 : 0.4;

      // Label dispatch: always-label set + hover.
      const showLabel = alwaysLabeledIds.has(n.id) || isHovered;

      // Draw node circle
      ctx.save();
      ctx.globalAlpha = nodeAlpha;
      ctx.beginPath();
      ctx.arc(node.x!, node.y!, r, 0, 2 * Math.PI);
      ctx.fillStyle = getNodeColor(n, colorBy);
      ctx.fill();

      // Outline: queried genes get a bold 3px ring; highlight/hover still shown.
      if (isQueried) {
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 3;
        ctx.stroke();
      } else if (isHighlighted || isHovered) {
        ctx.strokeStyle = isHighlighted ? "#ffffff" : "rgba(255,255,255,0.6)";
        ctx.lineWidth = isHighlighted ? 2 : 1;
        ctx.stroke();
      }
      ctx.restore();

      // Label (rendered at full alpha)
      if (showLabel) {
        const fontSize = Math.max(6, Math.min(10, 9 / globalScale));
        ctx.font = `600 ${fontSize}px Inter, system-ui, sans-serif`;
        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        ctx.fillStyle = "var(--color-foreground, #ffffff)";
        ctx.fillText(n.symbol, node.x!, node.y! + r + 1);
      }
    },
    [colorBy, highlightNodeId, hoveredNode, queriedIds, sharedCountOf, alwaysLabeledIds]
  );

  // Link rendering. Below-threshold edges dim to alpha 0.15 instead of being
  // filtered — that way users see the full structure and can adjust thresholds
  // without losing context.
  const paintLink = useCallback(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (link: any, ctx: CanvasRenderingContext2D) => {
      const edge = link as NetworkEdge;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const src = link.source as any;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const tgt = link.target as any;

      if (src.x == null || tgt.x == null) return;

      const { dim } = edgeVisible(edge, thresholds);

      ctx.beginPath();
      ctx.moveTo(src.x, src.y);
      ctx.lineTo(tgt.x, tgt.y);
      ctx.strokeStyle = getLayerColor(edge.layer);
      ctx.globalAlpha = dim ? 0.12 : 0.25 + edge.posterior * 0.35;
      ctx.lineWidth = dim ? 0.2 : 0.3 + edge.posterior * 0.8;
      ctx.stroke();
      ctx.globalAlpha = 1;
    },
    [thresholds]
  );

  // d3-force configuration
  useEffect(() => {
    const fg = fgRef.current;
    if (!fg) return;

    if (viewMode === "community") {
      // Disable force simulation for community map — use pre-computed FA2 coords
      fg.d3Force("charge", null);
      fg.d3Force("link", null);
      fg.d3Force("center", null);
    } else {
      // Neighborhood mode: stronger repulsion + longer links → more breathing
      // room between hubs and labels. Multi-gene queries loosen repulsion and
      // add an explicit center force so disconnected components don't drift
      // out of the viewport.
      fg.d3Force("charge")?.strength(multiGeneMode ? -90 : -180);
      fg.d3Force("link")?.distance(multiGeneMode ? 55 : 70);
      if (multiGeneMode) {
        fg.d3Force("center")?.strength?.(0.15);
      }
    }
  }, [viewMode, multiGeneMode]);

  // Reheat the simulation whenever the filtered node/link set changes so the
  // layout re-settles (prevents orphaned nodes from drifting away after filter
  // changes like F-stage / edge-type toggles).
  useEffect(() => {
    const fg = fgRef.current;
    if (!fg || viewMode === "community") return;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (fg as any).d3ReheatSimulation?.();
  }, [filteredNodes, filteredLinks, viewMode]);

  // Zoom to fit on data change
  useEffect(() => {
    const fg = fgRef.current;
    if (!fg || filteredNodes.length === 0) return;

    const timer = setTimeout(() => {
      fg.zoomToFit(400, 40);
    }, 500);

    return () => clearTimeout(timer);
  }, [filteredNodes.length]);

  // For community map mode, fix nodes to their FA2 positions
  const nodePositionFix = viewMode === "community"
    ? {
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        nodeVal: (n: any) => getNodeSize(n as NetworkNode) * 0.3,
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        d3VelocityDecay: 1, // Prevent movement
      }
    : {};

  return (
    <div ref={containerRef} className="relative h-full w-full overflow-hidden rounded-lg border border-border bg-card">
      {/* Empty state */}
      {filteredNodes.length === 0 && (
        <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">
          <div className="text-center">
            <svg className="mx-auto size-10 opacity-40" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1}>
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 3v2.25m6.364.386l-1.591 1.591M21 12h-2.25m-.386 6.364l-1.591-1.591M12 18.75V21m-4.773-4.227l-1.591 1.591M5.25 12H3m4.227-4.773L5.636 5.636M15.75 12a3.75 3.75 0 11-7.5 0 3.75 3.75 0 017.5 0z" />
            </svg>
            <p className="mt-2 text-sm">Search for a gene to explore its network</p>
          </div>
        </div>
      )}

      {filteredNodes.length > 0 && (
        <ForceGraph2D
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
          ref={fgRef as any}
          graphData={graphData}
          width={dimensions.w}
          height={dimensions.h}
          nodeCanvasObject={paintNode}
          nodePointerAreaPaint={(
            node: Record<string, unknown>,
            color: string,
            ctx: CanvasRenderingContext2D,
          ) => {
            const r = getNodeSize(node as unknown as NetworkNode);
            ctx.beginPath();
            ctx.arc(node.x as number, node.y as number, r + 2, 0, 2 * Math.PI);
            ctx.fillStyle = color;
            ctx.fill();
          }}
          linkCanvasObject={paintLink}
          onNodeClick={handleNodeClick}
          onNodeHover={handleNodeHover}
          cooldownTicks={viewMode === "community" ? 0 : 100}
          enableNodeDrag={viewMode !== "community"}
          backgroundColor="transparent"
          {...nodePositionFix}
        />
      )}

      {/* Stats overlay (bottom-left) */}
      {filteredNodes.length > 0 && (
        <div className="absolute bottom-3 left-3 rounded-md bg-card/80 px-2.5 py-1.5 text-[10px] text-muted-foreground backdrop-blur-sm">
          {filteredNodes.length.toLocaleString()} nodes &middot;{" "}
          {filteredLinks.length.toLocaleString()} edges
        </div>
      )}
    </div>
  );
}
