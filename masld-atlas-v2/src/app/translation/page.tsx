"use client";

import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import {
  ModalityMatrix,
  type ConvergenceRow,
} from "@/components/convergence/modality-matrix";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface GraphNode {
  id: string;
  type: "gene" | "drug" | "tf" | "pathway";
  label: string;
  // Gene-specific
  evidence_score?: number;
  is_deg?: boolean;
  is_coloc?: boolean;
  logFC?: number;
  // Drug-specific
  stage?: string;
  moa?: string;
  atlas_support?: string;
  // TF-specific
  n_targets?: number;
  activity_diff?: number;
  // Pathway-specific
  count?: number;
  // Computed layout
  x?: number;
  y?: number;
}

interface GraphEdge {
  source: string;
  target: string;
  type: "drug_target" | "regulon" | "pathway_member";
}

interface GraphData {
  nodes: GraphNode[];
  edges: GraphEdge[];
}

// ---------------------------------------------------------------------------
// Convergence stats
// ---------------------------------------------------------------------------

interface StatCard {
  value: string;
  label: string;
  sublabel: string;
}

const CONVERGENCE_STATS: StatCard[] = [
  {
    value: "7",
    label: "Modalities",
    sublabel: "Transcriptomic, genetic, epigenomic, spatial, single-cell, essentiality, drug",
  },
  {
    value: "2.6%",
    label: "DEG-COLOC overlap",
    sublabel: "Jaccard = 0.012, OR = 1.67",
  },
  {
    value: "123",
    label: "Pan-evidence genes",
    sublabel: "Active in all 7 sources",
  },
  {
    value: "18",
    label: "Triple-convergent",
    sublabel: "Transcriptomic + genetic + drug",
  },
];

// ---------------------------------------------------------------------------
// Layout constants
// ---------------------------------------------------------------------------

const SVG_W = 900;
const SVG_H = 640;
const PAD = 50;

// Regions for each node type
const GENE_CX = SVG_W / 2;
const GENE_CY = SVG_H / 2 - 10;
const GENE_RX = 230;
const GENE_RY = 200;

const DRUG_CX = 80;
const DRUG_CY = SVG_H / 2 - 10;

const TF_CX = SVG_W - 80;
const TF_CY = SVG_H / 2 - 10;

const PW_CX = SVG_W / 2;
const PW_CY = SVG_H - 50;

// ---------------------------------------------------------------------------
// Node colors
// ---------------------------------------------------------------------------

function geneColor(score: number): string {
  // Gradient: gray (0) -> blue (3) -> purple (6)
  if (score <= 0) return "hsl(220, 10%, 60%)";
  if (score === 1) return "hsl(220, 30%, 58%)";
  if (score === 2) return "hsl(225, 50%, 55%)";
  if (score === 3) return "hsl(230, 65%, 52%)";
  if (score === 4) return "hsl(250, 60%, 50%)";
  if (score === 5) return "hsl(265, 55%, 48%)";
  return "hsl(275, 60%, 45%)"; // 6+
}

const DRUG_COLOR = "hsl(152, 55%, 42%)";
const TF_COLOR = "hsl(30, 80%, 52%)";
const PW_COLOR = "hsl(220, 15%, 55%)";

const EDGE_COLORS: Record<string, string> = {
  drug_target: "hsl(152, 45%, 55%)",
  regulon: "hsl(30, 60%, 60%)",
  pathway_member: "hsl(220, 15%, 65%)",
};

// ---------------------------------------------------------------------------
// Layout algorithm
// ---------------------------------------------------------------------------

function layoutNodes(nodes: GraphNode[]): GraphNode[] {
  const positioned = nodes.map((n) => ({ ...n }));

  const genes = positioned.filter((n) => n.type === "gene");
  const drugs = positioned.filter((n) => n.type === "drug");
  const tfs = positioned.filter((n) => n.type === "tf");
  const pathways = positioned.filter((n) => n.type === "pathway");

  // Genes in an ellipse in the center
  genes.forEach((n, i) => {
    const angle = (2 * Math.PI * i) / genes.length - Math.PI / 2;
    n.x = GENE_CX + GENE_RX * Math.cos(angle);
    n.y = GENE_CY + GENE_RY * Math.sin(angle);
  });

  // Drugs on the left, spaced vertically
  const drugSpacing = Math.min(18, (SVG_H - 2 * PAD) / Math.max(drugs.length, 1));
  const drugStartY = DRUG_CY - ((drugs.length - 1) * drugSpacing) / 2;
  drugs.forEach((n, i) => {
    n.x = DRUG_CX;
    n.y = drugStartY + drugSpacing * i;
  });

  // TFs on the right, spaced vertically
  const tfSpacing = Math.min(22, (SVG_H - 2 * PAD) / Math.max(tfs.length, 1));
  const tfStartY = TF_CY - ((tfs.length - 1) * tfSpacing) / 2;
  tfs.forEach((n, i) => {
    n.x = TF_CX;
    n.y = tfStartY + tfSpacing * i;
  });

  // Pathways along the bottom
  const pwSpacing = Math.min(40, (SVG_W - 2 * PAD) / Math.max(pathways.length, 1));
  const pwStartX = PW_CX - ((pathways.length - 1) * pwSpacing) / 2;
  pathways.forEach((n, i) => {
    n.x = pwStartX + pwSpacing * i;
    n.y = PW_CY;
  });

  return positioned;
}

// ---------------------------------------------------------------------------
// Tooltip content builder
// ---------------------------------------------------------------------------

function tooltipContent(node: GraphNode): { title: string; lines: string[] } {
  const lines: string[] = [];
  if (node.type === "gene") {
    lines.push(`Evidence score: ${node.evidence_score ?? "N/A"}/7`);
    if (node.is_deg) lines.push("DEG: Yes");
    if (node.is_coloc) lines.push("COLOC: Yes");
    if (node.logFC != null) lines.push(`logFC: ${node.logFC >= 0 ? "+" : ""}${node.logFC.toFixed(3)}`);
  } else if (node.type === "drug") {
    if (node.moa) lines.push(`MoA: ${node.moa}`);
    if (node.stage) lines.push(node.stage);
    if (node.atlas_support) lines.push(`Atlas support: ${node.atlas_support}`);
  } else if (node.type === "tf") {
    if (node.n_targets != null) lines.push(`Targets in graph: ${node.n_targets}`);
    if (node.activity_diff != null) lines.push(`Activity diff: ${node.activity_diff >= 0 ? "+" : ""}${node.activity_diff.toFixed(4)}`);
  } else if (node.type === "pathway") {
    if (node.count != null) lines.push(`Genes in pathway: ${node.count}`);
  }
  return { title: node.label, lines };
}

// ---------------------------------------------------------------------------
// SVG node shapes
// ---------------------------------------------------------------------------

function GeneCircle({
  node,
  r,
  isHighlighted,
  isDimmed,
  onMouseEnter,
  onMouseLeave,
  onClick,
}: {
  node: GraphNode;
  r: number;
  isHighlighted: boolean;
  isDimmed: boolean;
  onMouseEnter: () => void;
  onMouseLeave: () => void;
  onClick: () => void;
}) {
  return (
    <circle
      cx={node.x}
      cy={node.y}
      r={isHighlighted ? r + 2 : r}
      fill={geneColor(node.evidence_score ?? 0)}
      opacity={isDimmed ? 0.15 : isHighlighted ? 1 : 0.8}
      stroke={isHighlighted ? "hsl(0, 0%, 100%)" : node.is_coloc ? "hsl(45, 90%, 60%)" : "none"}
      strokeWidth={isHighlighted ? 2 : node.is_coloc ? 1.5 : 0}
      style={{ cursor: "pointer", transition: "opacity 0.15s, r 0.15s" }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
      onClick={onClick}
    />
  );
}

function DrugDiamond({
  node,
  size,
  isHighlighted,
  isDimmed,
  onMouseEnter,
  onMouseLeave,
}: {
  node: GraphNode;
  size: number;
  isHighlighted: boolean;
  isDimmed: boolean;
  onMouseEnter: () => void;
  onMouseLeave: () => void;
}) {
  const s = isHighlighted ? size + 2 : size;
  const x = node.x ?? 0;
  const y = node.y ?? 0;
  const points = `${x},${y - s} ${x + s},${y} ${x},${y + s} ${x - s},${y}`;
  return (
    <polygon
      points={points}
      fill={DRUG_COLOR}
      opacity={isDimmed ? 0.15 : isHighlighted ? 1 : 0.8}
      stroke={isHighlighted ? "hsl(0, 0%, 100%)" : "none"}
      strokeWidth={isHighlighted ? 2 : 0}
      style={{ cursor: "pointer", transition: "opacity 0.15s" }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
    />
  );
}

function TfTriangle({
  node,
  size,
  isHighlighted,
  isDimmed,
  onMouseEnter,
  onMouseLeave,
}: {
  node: GraphNode;
  size: number;
  isHighlighted: boolean;
  isDimmed: boolean;
  onMouseEnter: () => void;
  onMouseLeave: () => void;
}) {
  const s = isHighlighted ? size + 2 : size;
  const x = node.x ?? 0;
  const y = node.y ?? 0;
  const points = `${x},${y - s} ${x + s * 0.87},${y + s * 0.5} ${x - s * 0.87},${y + s * 0.5}`;
  return (
    <polygon
      points={points}
      fill={TF_COLOR}
      opacity={isDimmed ? 0.15 : isHighlighted ? 1 : 0.8}
      stroke={isHighlighted ? "hsl(0, 0%, 100%)" : "none"}
      strokeWidth={isHighlighted ? 2 : 0}
      style={{ cursor: "pointer", transition: "opacity 0.15s" }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
    />
  );
}

function PathwayRect({
  node,
  w,
  h,
  isHighlighted,
  isDimmed,
  onMouseEnter,
  onMouseLeave,
}: {
  node: GraphNode;
  w: number;
  h: number;
  isHighlighted: boolean;
  isDimmed: boolean;
  onMouseEnter: () => void;
  onMouseLeave: () => void;
}) {
  const x = (node.x ?? 0) - w / 2;
  const y = (node.y ?? 0) - h / 2;
  return (
    <rect
      x={x}
      y={y}
      width={w}
      height={h}
      rx={4}
      fill={PW_COLOR}
      opacity={isDimmed ? 0.15 : isHighlighted ? 1 : 0.7}
      stroke={isHighlighted ? "hsl(0, 0%, 100%)" : "none"}
      strokeWidth={isHighlighted ? 2 : 0}
      style={{ cursor: "pointer", transition: "opacity 0.15s" }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
    />
  );
}

// ---------------------------------------------------------------------------
// Knowledge Graph component
// ---------------------------------------------------------------------------

type NodeType = "gene" | "drug" | "tf" | "pathway";

interface KnowledgeGraphProps {
  data: GraphData;
}

function KnowledgeGraph({ data }: KnowledgeGraphProps) {
  const router = useRouter();
  const [hoveredId, setHoveredId] = useState<string | null>(null);
  const [visibleTypes, setVisibleTypes] = useState<Record<NodeType, boolean>>({
    gene: true,
    drug: true,
    tf: true,
    pathway: true,
  });

  // Layout nodes once
  const positioned = useMemo(() => layoutNodes(data.nodes), [data.nodes]);

  // Build adjacency map for hover highlighting
  const adjacency = useMemo(() => {
    const map = new Map<string, Set<string>>();
    for (const edge of data.edges) {
      if (!map.has(edge.source)) map.set(edge.source, new Set());
      if (!map.has(edge.target)) map.set(edge.target, new Set());
      map.get(edge.source)!.add(edge.target);
      map.get(edge.target)!.add(edge.source);
    }
    return map;
  }, [data.edges]);

  // Build node lookup
  const nodeMap = useMemo(() => {
    const m = new Map<string, GraphNode>();
    for (const n of positioned) m.set(n.id, n);
    return m;
  }, [positioned]);

  // Determine which nodes are highlighted / dimmed
  const connectedIds = useMemo(() => {
    if (!hoveredId) return null;
    const set = new Set<string>();
    set.add(hoveredId);
    const neighbors = adjacency.get(hoveredId);
    if (neighbors) neighbors.forEach((nid) => set.add(nid));
    return set;
  }, [hoveredId, adjacency]);

  // Filter visible nodes and edges
  const visibleNodes = useMemo(
    () => positioned.filter((n) => visibleTypes[n.type]),
    [positioned, visibleTypes]
  );
  const visibleNodeIds = useMemo(
    () => new Set(visibleNodes.map((n) => n.id)),
    [visibleNodes]
  );
  const visibleEdges = useMemo(
    () => data.edges.filter((e) => visibleNodeIds.has(e.source) && visibleNodeIds.has(e.target)),
    [data.edges, visibleNodeIds]
  );

  function isHighlighted(id: string): boolean {
    if (!connectedIds) return false;
    return connectedIds.has(id);
  }
  function isDimmed(id: string): boolean {
    if (!connectedIds) return false;
    return !connectedIds.has(id);
  }

  // Tooltip state
  const hoveredNode = hoveredId ? nodeMap.get(hoveredId) ?? null : null;
  const tip = hoveredNode ? tooltipContent(hoveredNode) : null;

  function toggleType(t: NodeType) {
    setVisibleTypes((prev) => ({ ...prev, [t]: !prev[t] }));
  }

  // Compute tooltip position
  const tipX = hoveredNode?.x != null ? Math.min(hoveredNode.x + 14, SVG_W - 180) : 0;
  const tipY = hoveredNode?.y != null ? Math.max(hoveredNode.y - 20, 10) : 0;

  return (
    <div>
      {/* Filter toggles */}
      <div className="mb-3 flex flex-wrap items-center gap-2">
        <span className="text-xs font-medium text-muted-foreground">Show:</span>
        {([
          { type: "gene" as NodeType, label: "Genes", color: "hsl(275, 60%, 45%)", shape: "circle" },
          { type: "drug" as NodeType, label: "Drugs", color: DRUG_COLOR, shape: "diamond" },
          { type: "tf" as NodeType, label: "TFs", color: TF_COLOR, shape: "triangle" },
          { type: "pathway" as NodeType, label: "Pathways", color: PW_COLOR, shape: "rect" },
        ]).map((item) => (
          <button
            key={item.type}
            onClick={() => toggleType(item.type)}
            className={`flex items-center gap-1.5 rounded-md border px-2.5 py-1 text-xs font-medium transition-colors ${
              visibleTypes[item.type]
                ? "border-border bg-card text-foreground"
                : "border-border/50 bg-muted/30 text-muted-foreground line-through"
            }`}
          >
            <span
              className="inline-block h-2.5 w-2.5 rounded-sm"
              style={{ background: item.color, opacity: visibleTypes[item.type] ? 1 : 0.3 }}
            />
            {item.label} ({positioned.filter((n) => n.type === item.type).length})
          </button>
        ))}
      </div>

      {/* SVG graph */}
      <div className="overflow-x-auto rounded-lg border border-border bg-card shadow-sm">
        <svg
          width={SVG_W}
          height={SVG_H}
          viewBox={`0 0 ${SVG_W} ${SVG_H}`}
          className="max-w-full"
          style={{ fontFamily: "inherit" }}
        >
          {/* Background */}
          <rect width={SVG_W} height={SVG_H} className="fill-muted/20" rx={8} />

          {/* Region labels */}
          <text x={DRUG_CX} y={30} textAnchor="middle" fontSize={10} fontWeight={600} fill="hsl(var(--muted-foreground))" opacity={0.7}>
            DRUGS
          </text>
          <text x={TF_CX} y={30} textAnchor="middle" fontSize={10} fontWeight={600} fill="hsl(var(--muted-foreground))" opacity={0.7}>
            TRANSCRIPTION FACTORS
          </text>
          <text x={PW_CX} y={SVG_H - 12} textAnchor="middle" fontSize={10} fontWeight={600} fill="hsl(var(--muted-foreground))" opacity={0.7}>
            PATHWAYS
          </text>
          <text x={GENE_CX} y={GENE_CY - GENE_RY - 18} textAnchor="middle" fontSize={10} fontWeight={600} fill="hsl(var(--muted-foreground))" opacity={0.7}>
            GENES ({positioned.filter((n) => n.type === "gene").length})
          </text>

          {/* Edges */}
          {visibleEdges.map((edge, i) => {
            const src = nodeMap.get(edge.source);
            const tgt = nodeMap.get(edge.target);
            if (!src?.x || !tgt?.x || src.y == null || tgt.y == null) return null;
            const highlighted = hoveredId != null && (hoveredId === edge.source || hoveredId === edge.target);
            const dimmed = hoveredId != null && !highlighted;
            return (
              <line
                key={`e-${i}`}
                x1={src.x}
                y1={src.y}
                x2={tgt.x}
                y2={tgt.y}
                stroke={EDGE_COLORS[edge.type] ?? "hsl(var(--border))"}
                strokeWidth={highlighted ? 1.8 : 0.6}
                opacity={dimmed ? 0.05 : highlighted ? 0.9 : 0.2}
                style={{ transition: "opacity 0.15s, stroke-width 0.15s" }}
              />
            );
          })}

          {/* Nodes */}
          {visibleNodes.map((node) => {
            const hl = isHighlighted(node.id);
            const dm = isDimmed(node.id);
            const enterHandler = () => setHoveredId(node.id);
            const leaveHandler = () => setHoveredId(null);

            if (node.type === "gene") {
              return (
                <GeneCircle
                  key={node.id}
                  node={node}
                  r={node.is_coloc ? 5.5 : 4}
                  isHighlighted={hl}
                  isDimmed={dm}
                  onMouseEnter={enterHandler}
                  onMouseLeave={leaveHandler}
                  onClick={() => router.push(`/gene/${encodeURIComponent(node.label)}`)}
                />
              );
            }
            if (node.type === "drug") {
              return (
                <DrugDiamond
                  key={node.id}
                  node={node}
                  size={6}
                  isHighlighted={hl}
                  isDimmed={dm}
                  onMouseEnter={enterHandler}
                  onMouseLeave={leaveHandler}
                />
              );
            }
            if (node.type === "tf") {
              return (
                <TfTriangle
                  key={node.id}
                  node={node}
                  size={7}
                  isHighlighted={hl}
                  isDimmed={dm}
                  onMouseEnter={enterHandler}
                  onMouseLeave={leaveHandler}
                />
              );
            }
            if (node.type === "pathway") {
              return (
                <PathwayRect
                  key={node.id}
                  node={node}
                  w={12}
                  h={10}
                  isHighlighted={hl}
                  isDimmed={dm}
                  onMouseEnter={enterHandler}
                  onMouseLeave={leaveHandler}
                />
              );
            }
            return null;
          })}

          {/* Node labels (only for highlighted or always-labeled nodes) */}
          {visibleNodes.map((node) => {
            if (!node.x || node.y == null) return null;
            const hl = isHighlighted(node.id);
            // Always show drug/TF labels; gene/pathway labels only on hover
            const showLabel = node.type === "drug" || node.type === "tf" || hl;
            if (!showLabel) return null;

            const labelX = node.type === "drug" ? node.x - 12 : node.type === "tf" ? node.x + 12 : node.x;
            const labelY = node.type === "drug" || node.type === "tf" ? node.y + 3 : node.y - 10;
            const anchor = node.type === "drug" ? "end" : node.type === "tf" ? "start" : "middle";

            return (
              <text
                key={`lbl-${node.id}`}
                x={labelX}
                y={labelY}
                textAnchor={anchor}
                fontSize={node.type === "gene" ? 9 : 8}
                fontWeight={hl ? 700 : 500}
                fill="hsl(var(--foreground))"
                opacity={hl ? 1 : 0.7}
                style={{ pointerEvents: "none", transition: "opacity 0.15s" }}
              >
                {node.label.length > 20 ? node.label.slice(0, 18) + "..." : node.label}
              </text>
            );
          })}

          {/* Tooltip */}
          {tip && hoveredNode?.x != null && hoveredNode?.y != null && (
            <g style={{ pointerEvents: "none" }}>
              <rect
                x={tipX}
                y={tipY}
                width={170}
                height={20 + tip.lines.length * 14}
                rx={4}
                fill="hsl(var(--popover))"
                stroke="hsl(var(--border))"
                strokeWidth={1}
              />
              <text x={tipX + 8} y={tipY + 15} fontSize={11} fontWeight={700} fill="hsl(var(--foreground))">
                {tip.title}
              </text>
              {tip.lines.map((line, i) => (
                <text
                  key={i}
                  x={tipX + 8}
                  y={tipY + 29 + i * 14}
                  fontSize={10}
                  fill="hsl(var(--muted-foreground))"
                >
                  {line}
                </text>
              ))}
            </g>
          )}
        </svg>
      </div>

      {/* Legend row */}
      <div className="mt-2 flex flex-wrap gap-4 text-xs text-muted-foreground">
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: "hsl(275, 60%, 45%)" }} />
          Gene (score 6-7)
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: "hsl(230, 65%, 52%)" }} />
          Gene (score 3)
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: "hsl(220, 10%, 60%)" }} />
          Gene (score 0)
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full" style={{ background: "hsl(45, 90%, 60%)", border: "1.5px solid hsl(45, 90%, 60%)" }} />
          COLOC ring
        </span>
        <span className="ml-auto italic">Click a gene node to view profile</span>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Bayesian integration results
// ---------------------------------------------------------------------------

interface BayesianGene {
  rank: number;
  symbol: string;
  score: number;
  percentile: number;
  is_deg: boolean;
  is_coloc: boolean;
  is_druggable: boolean;
}

interface BayesianRanking {
  genes: BayesianGene[];
  metadata: { n_genes: number; score_source: string; notes?: string };
}

const HIGHLIGHT_GENES = ["THRB", "NR1H4", "PPARA"];

type BayesSortKey = "rank" | "symbol" | "score" | "percentile";

function BayesianRankingTable({ data }: { data: BayesianRanking }) {
  const [query, setQuery] = useState("");
  const [sortKey, setSortKey] = useState<BayesSortKey>("rank");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");
  const [page, setPage] = useState(0);
  const PAGE_SIZE = 100;

  const filtered = useMemo(() => {
    const q = query.trim().toUpperCase();
    let rows = data.genes;
    if (q) {
      rows = rows.filter((g) => g.symbol.toUpperCase().includes(q));
    }
    const sorted = [...rows];
    sorted.sort((a, b) => {
      const av = a[sortKey];
      const bv = b[sortKey];
      let cmp = 0;
      if (typeof av === "number" && typeof bv === "number") cmp = av - bv;
      else cmp = String(av).localeCompare(String(bv));
      return sortDir === "asc" ? cmp : -cmp;
    });
    return sorted;
  }, [data, query, sortKey, sortDir]);

  const nPages = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE));
  const pageSafe = Math.min(page, nPages - 1);
  const pageRows = filtered.slice(
    pageSafe * PAGE_SIZE,
    (pageSafe + 1) * PAGE_SIZE
  );

  const toggleSort = (k: BayesSortKey) => {
    if (k === sortKey) setSortDir(sortDir === "asc" ? "desc" : "asc");
    else {
      setSortKey(k);
      setSortDir(k === "rank" ? "asc" : "desc");
    }
    setPage(0);
  };

  const pinned = useMemo(
    () =>
      HIGHLIGHT_GENES
        .map((sym) => data.genes.find((g) => g.symbol === sym))
        .filter(Boolean) as BayesianGene[],
    [data]
  );

  const headers: { key: BayesSortKey; label: string }[] = [
    { key: "rank", label: "Rank" },
    { key: "symbol", label: "Gene" },
    { key: "score", label: "Score" },
    { key: "percentile", label: "Percentile" },
  ];

  return (
    <div className="flex flex-col gap-3">
      {/* Pinned callouts */}
      {pinned.length > 0 && (
        <div className="rounded-md border border-primary/30 bg-primary/5 px-4 py-2.5 text-sm">
          <span className="mr-2 text-xs font-medium text-muted-foreground">
            Highlighted:
          </span>
          {pinned.map((g, i) => (
            <span key={g.symbol} className="font-mono text-xs">
              {i > 0 && <span className="mx-1.5 text-muted-foreground">·</span>}
              <Link
                href={`/gene/${encodeURIComponent(g.symbol)}`}
                className="font-semibold text-primary hover:underline"
              >
                {g.symbol}
              </Link>{" "}
              <span className="text-muted-foreground">
                (rank {g.rank.toLocaleString()}, top {g.percentile.toFixed(1)}%)
              </span>
            </span>
          ))}
        </div>
      )}

      {/* Controls */}
      <div className="flex flex-wrap items-center gap-3 text-sm">
        <input
          type="search"
          placeholder="Search gene…"
          value={query}
          onChange={(e) => {
            setQuery(e.target.value);
            setPage(0);
          }}
          className="h-8 w-48 rounded-md border border-border bg-background px-2 text-sm"
        />
        <span className="text-xs text-muted-foreground">
          {filtered.length.toLocaleString()} of{" "}
          {data.metadata.n_genes.toLocaleString()} genes
        </span>
        <div className="ml-auto flex items-center gap-2">
          <button
            type="button"
            onClick={() => setPage((p) => Math.max(0, p - 1))}
            disabled={pageSafe === 0}
            className="rounded-md border border-border px-2 py-1 text-xs hover:bg-muted disabled:opacity-40"
          >
            ← Prev
          </button>
          <span className="font-mono text-xs text-muted-foreground">
            {pageSafe + 1} / {nPages}
          </span>
          <button
            type="button"
            onClick={() => setPage((p) => Math.min(nPages - 1, p + 1))}
            disabled={pageSafe >= nPages - 1}
            className="rounded-md border border-border px-2 py-1 text-xs hover:bg-muted disabled:opacity-40"
          >
            Next →
          </button>
        </div>
      </div>

      {/* Table */}
      <div className="overflow-x-auto rounded-lg border border-border">
        <table className="w-full text-sm">
          <thead className="sticky top-0 bg-muted/80 backdrop-blur">
            <tr className="border-b border-border">
              {headers.map((h) => (
                <th
                  key={h.key}
                  onClick={() => toggleSort(h.key)}
                  className="cursor-pointer px-4 py-2.5 text-left text-xs font-medium text-muted-foreground hover:text-foreground"
                >
                  {h.label}
                  {sortKey === h.key && (
                    <span className="ml-1">
                      {sortDir === "asc" ? "▲" : "▼"}
                    </span>
                  )}
                </th>
              ))}
              <th className="px-4 py-2.5 text-left text-xs font-medium text-muted-foreground">
                DEG
              </th>
              <th className="px-4 py-2.5 text-left text-xs font-medium text-muted-foreground">
                COLOC
              </th>
              <th className="px-4 py-2.5 text-left text-xs font-medium text-muted-foreground">
                Druggable
              </th>
            </tr>
          </thead>
          <tbody>
            {pageRows.map((g) => (
              <tr
                key={g.symbol}
                className="border-b border-border/50 transition-colors hover:bg-muted/30"
              >
                <td className="px-4 py-2 font-mono text-xs">
                  #{g.rank.toLocaleString()}
                </td>
                <td className="px-4 py-2">
                  <Link
                    href={`/gene/${encodeURIComponent(g.symbol)}`}
                    className="font-mono text-xs font-semibold text-primary hover:underline"
                  >
                    {g.symbol}
                  </Link>
                </td>
                <td className="px-4 py-2 font-mono text-xs">
                  {g.score.toFixed(4)}
                </td>
                <td className="px-4 py-2 font-mono text-xs text-muted-foreground">
                  {g.percentile.toFixed(2)}%
                </td>
                <td className="px-4 py-2 text-xs">
                  {g.is_deg ? "✓" : ""}
                </td>
                <td className="px-4 py-2 text-xs">
                  {g.is_coloc ? "✓" : ""}
                </td>
                <td className="px-4 py-2 text-xs">
                  {g.is_druggable ? "✓" : ""}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p className="text-[11px] text-muted-foreground">
        Score source: {data.metadata.score_source}
      </p>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

type ConvergenceTab = "matrix" | "graph" | "bayesian";

export default function TranslationPage() {
  const [graphData, setGraphData] = useState<GraphData | null>(null);
  const [loading, setLoading] = useState(true);
  const [convergence, setConvergence] = useState<ConvergenceRow[] | null>(null);
  const [bayesian, setBayesian] = useState<BayesianRanking | null>(null);
  const [tab, setTab] = useState<ConvergenceTab>("matrix");

  useEffect(() => {
    fetch("/data/knowledge_graph.json")
      .then((res) => res.json())
      .then((data: GraphData) => {
        setGraphData(data);
        setLoading(false);
      })
      .catch(() => setLoading(false));
    fetch("/data/convergence_matrix.json")
      .then((res) => res.json())
      .then((data: ConvergenceRow[]) => setConvergence(data))
      .catch(() => setConvergence([]));
    fetch("/data/bayesian_ranking.json")
      .then((res) => res.json())
      .then((data: BayesianRanking) => setBayesian(data))
      .catch(() =>
        setBayesian({
          genes: [],
          metadata: { n_genes: 0, score_source: "unavailable" },
        })
      );
  }, []);

  const TABS: { key: ConvergenceTab; label: string }[] = [
    { key: "matrix", label: "Convergence Matrix" },
    { key: "graph", label: "Knowledge Graph" },
    { key: "bayesian", label: "Bayesian Ranking" },
  ];

  return (
    <div className="mx-auto max-w-5xl px-6 py-10">
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Convergence &amp; Translation
        </h1>
        <p className="mt-2 text-muted-foreground">
          Multi-evidence integration reveals that each source captures
          non-overlapping biology — convergence, not any single analysis,
          most reliably identifies therapeutic targets.
        </p>
      </div>

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Convergence Overview                                      */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Convergence Overview
        </h2>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
          {CONVERGENCE_STATS.map((card) => (
            <div
              key={card.label}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm"
            >
              <p className="font-mono text-2xl font-bold text-primary">
                {card.value}
              </p>
              <p className="mt-0.5 text-sm font-semibold">{card.label}</p>
              <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                {card.sublabel}
              </p>
            </div>
          ))}
        </div>
        <p className="mt-4 text-sm text-muted-foreground">
          Only 2.6% of DEGs have colocalization support, demonstrating that
          transcriptomic and genetic evidence capture fundamentally different
          biology. Multi-evidence convergence across 7 independent sources
          provides the most reliable path to therapeutic target identification.
        </p>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Tabbed convergence views                                             */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <div className="mb-3 flex items-center gap-1 border-b border-border">
          {TABS.map((t) => (
            <button
              key={t.key}
              type="button"
              onClick={() => setTab(t.key)}
              className={
                "rounded-t-md border-b-2 px-4 py-2 text-sm font-medium transition-colors " +
                (tab === t.key
                  ? "border-primary text-primary"
                  : "border-transparent text-muted-foreground hover:text-foreground")
              }
            >
              {t.label}
            </button>
          ))}
        </div>

        {tab === "matrix" && (
          <div className="pt-4">
            {convergence === null ? (
              <div className="flex h-[500px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
                Loading convergence matrix&hellip;
              </div>
            ) : convergence.length === 0 ? (
              <div className="flex h-[500px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
                Failed to load convergence_matrix.json
              </div>
            ) : (
              <ModalityMatrix data={convergence} />
            )}
          </div>
        )}

        {tab === "graph" && (
          <div className="pt-4">
            <div className="mb-3 flex items-center justify-between">
              <span className="text-xs text-muted-foreground">
                Knowledge graph: gene × drug × TF × pathway
              </span>
              {graphData && (
                <span className="text-xs text-muted-foreground">
                  {graphData.nodes.length} nodes, {graphData.edges.length} edges
                </span>
              )}
            </div>
            {loading ? (
              <div className="flex h-[640px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
                Loading knowledge graph&hellip;
              </div>
            ) : graphData ? (
              <KnowledgeGraph data={graphData} />
            ) : (
              <div className="flex h-[640px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
                Failed to load knowledge graph data.
              </div>
            )}
            <p className="mt-2 text-xs text-muted-foreground">
              Genes are arranged in a central ellipse colored by evidence score
              (gray = low, purple = high). Gold ring indicates COLOC support.
              Drugs (left), transcription factors (right), and pathways
              (bottom) are connected by their respective relationships. Hover
              to highlight connections; click any gene to view its profile.
            </p>
          </div>
        )}

        {tab === "bayesian" && (
          <div className="pt-4">
            <p className="mb-4 text-sm text-muted-foreground">
              Cumulative enrichment, Bayesian posterior integration, and
              network propagation (Scripts 46a/46b/46c) combine all 7
              modalities into a unified gene ranking. Rank 1 = most evidence.
            </p>
            {bayesian === null ? (
              <div className="flex h-[400px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
                Loading Bayesian ranking&hellip;
              </div>
            ) : (
              <BayesianRankingTable data={bayesian} />
            )}
          </div>
        )}
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: Modality Independence                                      */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Modality Independence
        </h2>
        <div className="rounded-lg border border-border bg-card px-5 py-4 shadow-sm">
          <p className="text-sm leading-relaxed text-foreground">
            196 genes are both differentially expressed and have colocalization
            support (OR&nbsp;=&nbsp;1.67, P&nbsp;=&nbsp;5.5&times;10
            <sup>-6</sup>). While statistically significant, the Jaccard
            index of 0.012 confirms that the two evidence layers are nearly
            non-overlapping. This low overlap is not a weakness — it means
            each modality captures biology invisible to the other.
          </p>
          <div className="mt-3 flex flex-wrap gap-2">
            <Badge variant="outline" className="text-[10px]">
              5,484 DEGs (transcriptomic)
            </Badge>
            <Badge variant="outline" className="text-[10px]">
              348 COLOC genes (genetic)
            </Badge>
            <Badge variant="outline" className="text-[10px]">
              196 overlap (2.6%)
            </Badge>
            <Badge variant="outline" className="text-[10px]">
              OR = 1.67
            </Badge>
          </div>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/causal"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          &larr; Causal Architecture
        </Link>
        <Link
          href="/drugs"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Drug Pipeline &rarr;
        </Link>
        <Link
          href="/explore"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Browse DEGs &rarr;
        </Link>
      </div>
    </div>
  );
}
