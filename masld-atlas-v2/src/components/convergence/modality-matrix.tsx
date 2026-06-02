"use client";

import { useMemo, useState } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export interface ConvergenceRow {
  gene: string;
  modalities: number[]; // length 8
  count: number;
  is_deg: boolean;
  is_coloc: boolean;
  is_druggable: boolean;
}

interface ModalityMatrixProps {
  data: ConvergenceRow[];
  onRowClick?: (gene: string) => void;
}

// Modality headers — order matches `modalities[]`
const MODALITY_HEADERS = [
  { key: "m1", label: "M1", title: "Human bulk DE (Integrated logFC)" },
  { key: "m2", label: "M2", title: "Mouse conserved (meta-analysis padj)" },
  { key: "m3", label: "M3", title: "Genetic causal (SuSiE-COLOC PP4)" },
  { key: "m4", label: "M4", title: "Essentiality (DepMap CHRONOS)" },
  { key: "m5", label: "M5", title: "Epigenomic (GWAS-ATAC regulatory)" },
  { key: "m6", label: "M6", title: "Spatial (Visium Moran's I / SVG)" },
  { key: "m7", label: "M7", title: "Single-cell (cross-cell-type DE / LIANA)" },
  { key: "m8", label: "M8", title: "Proteomics (best protein logFC)" },
];

const ROW_HEIGHT = 18;
const ROW_LABEL_WIDTH = 92;
const COL_WIDTH = 44;
const HEADER_HEIGHT = 48;
const MATRIX_HEIGHT = 580;
const DEFAULT_TOP_N = 500;

// Single-hue orange scale: 0 → white, 1 → dark orange
function cellColor(v: number): string {
  if (v <= 0) return "hsl(0 0% 97%)";
  const clamped = Math.min(Math.max(v, 0), 1);
  // lightness 95% → 45%
  const l = 95 - clamped * 50;
  return `hsl(24 90% ${l.toFixed(1)}%)`;
}

type SortMode = "count" | "alpha" | "m1" | "m3";

export function ModalityMatrix({ data, onRowClick }: ModalityMatrixProps) {
  const [query, setQuery] = useState("");
  const [minCount, setMinCount] = useState(0);
  const [requiredModality, setRequiredModality] = useState<number | null>(null);
  const [sortMode, setSortMode] = useState<SortMode>("count");
  const [showAll, setShowAll] = useState(false);

  const filtered = useMemo(() => {
    const q = query.trim().toUpperCase();
    let rows = data;
    if (q) {
      rows = rows.filter((r) => r.gene.toUpperCase().includes(q));
    }
    if (minCount > 0) {
      rows = rows.filter((r) => r.count >= minCount);
    }
    if (requiredModality !== null) {
      rows = rows.filter((r) => (r.modalities[requiredModality] ?? 0) > 0.1);
    }
    // Sort
    const sorted = [...rows];
    switch (sortMode) {
      case "alpha":
        sorted.sort((a, b) => a.gene.localeCompare(b.gene));
        break;
      case "m1":
        sorted.sort((a, b) => (b.modalities[0] ?? 0) - (a.modalities[0] ?? 0));
        break;
      case "m3":
        sorted.sort((a, b) => (b.modalities[2] ?? 0) - (a.modalities[2] ?? 0));
        break;
      case "count":
      default:
        sorted.sort((a, b) => b.count - a.count || a.gene.localeCompare(b.gene));
    }
    return sorted;
  }, [data, query, minCount, requiredModality, sortMode]);

  const totalGenes = data.length;
  const displayedRows = showAll ? filtered : filtered.slice(0, DEFAULT_TOP_N);
  const maxCount = useMemo(
    () => data.reduce((m, r) => (r.count > m ? r.count : m), 0),
    [data]
  );
  const countHistogram = useMemo(() => {
    const h = new Array(8).fill(0);
    data.forEach((r) => {
      h[Math.min(r.count, 7)] += 1;
    });
    return h;
  }, [data]);

  const svgHeight = HEADER_HEIGHT + displayedRows.length * ROW_HEIGHT + 4;

  return (
    <div className="flex flex-col gap-3">
      {/* Summary */}
      <div className="rounded-md border border-border bg-card px-4 py-3 text-sm">
        <p>
          Showing{" "}
          <span className="font-mono font-semibold">
            {displayedRows.length.toLocaleString()}
          </span>{" "}
          of{" "}
          <span className="font-mono font-semibold">
            {totalGenes.toLocaleString()}
          </span>{" "}
          genes. Highest convergence:{" "}
          <span className="font-mono font-semibold">{maxCount}</span>/7
          modalities.
        </p>
        <p className="mt-1 text-xs text-muted-foreground">
          By modality count:{" "}
          {countHistogram
            .map((n, i) => (n > 0 ? `${i}/7=${n.toLocaleString()}` : null))
            .filter(Boolean)
            .join(" · ")}
          .
        </p>
      </div>

      {/* Controls */}
      <div className="flex flex-wrap items-center gap-3 text-sm">
        <input
          type="search"
          placeholder="Search gene…"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          className="h-8 w-48 rounded-md border border-border bg-background px-2 text-sm"
        />
        <label className="flex items-center gap-2">
          <span className="text-xs text-muted-foreground">Min modalities:</span>
          <input
            type="range"
            min={0}
            max={7}
            step={1}
            value={minCount}
            onChange={(e) => setMinCount(parseInt(e.target.value))}
            className="w-32"
          />
          <span className="font-mono text-xs">{minCount}</span>
        </label>
        <label className="flex items-center gap-2">
          <span className="text-xs text-muted-foreground">Sort:</span>
          <select
            value={sortMode}
            onChange={(e) => setSortMode(e.target.value as SortMode)}
            className="h-8 rounded-md border border-border bg-background px-2 text-xs"
          >
            <option value="count">Modality count</option>
            <option value="alpha">Gene (A–Z)</option>
            <option value="m1">M1 Human DE</option>
            <option value="m3">M3 Genetic (COLOC)</option>
          </select>
        </label>
        {requiredModality !== null && (
          <Badge variant="outline" className="text-[10px]">
            M{requiredModality + 1} required
            <button
              type="button"
              onClick={() => setRequiredModality(null)}
              className="ml-1 text-muted-foreground hover:text-foreground"
            >
              ×
            </button>
          </Badge>
        )}
        <div className="ml-auto flex items-center gap-2 text-xs text-muted-foreground">
          {!showAll && filtered.length > DEFAULT_TOP_N && (
            <button
              type="button"
              onClick={() => setShowAll(true)}
              className="rounded-md border border-border px-2 py-1 hover:bg-muted"
            >
              Show all {filtered.length.toLocaleString()}
            </button>
          )}
          {showAll && (
            <button
              type="button"
              onClick={() => setShowAll(false)}
              className="rounded-md border border-border px-2 py-1 hover:bg-muted"
            >
              Top {DEFAULT_TOP_N}
            </button>
          )}
        </div>
      </div>

      {/* Matrix (virtualized-lite: scrollable container) */}
      <div
        className="relative overflow-auto rounded-md border border-border bg-background"
        style={{ maxHeight: MATRIX_HEIGHT }}
      >
        <svg
          width={ROW_LABEL_WIDTH + MODALITY_HEADERS.length * COL_WIDTH + 40}
          height={svgHeight}
          style={{ display: "block" }}
        >
          {/* Header row */}
          <g>
            <rect
              x={0}
              y={0}
              width={ROW_LABEL_WIDTH + MODALITY_HEADERS.length * COL_WIDTH + 40}
              height={HEADER_HEIGHT}
              fill="hsl(0 0% 99%)"
              stroke="hsl(0 0% 90%)"
            />
            {MODALITY_HEADERS.map((m, i) => {
              const active = requiredModality === i;
              return (
                <g
                  key={m.key}
                  transform={`translate(${ROW_LABEL_WIDTH + i * COL_WIDTH},0)`}
                  onClick={() =>
                    setRequiredModality(active ? null : i)
                  }
                  style={{ cursor: "pointer" }}
                >
                  <title>{m.title} (click to toggle filter)</title>
                  <rect
                    x={1}
                    y={1}
                    width={COL_WIDTH - 2}
                    height={HEADER_HEIGHT - 2}
                    fill={active ? "hsl(24 90% 85%)" : "transparent"}
                  />
                  <text
                    x={COL_WIDTH / 2}
                    y={HEADER_HEIGHT / 2 + 4}
                    textAnchor="middle"
                    fontSize={11}
                    fontFamily="monospace"
                    fontWeight={600}
                    fill="hsl(0 0% 25%)"
                  >
                    {m.label}
                  </text>
                </g>
              );
            })}
            <text
              x={ROW_LABEL_WIDTH + MODALITY_HEADERS.length * COL_WIDTH + 6}
              y={HEADER_HEIGHT / 2 + 4}
              fontSize={10}
              fill="hsl(0 0% 45%)"
            >
              #
            </text>
          </g>
          {/* Rows */}
          {displayedRows.map((r, idx) => {
            const y = HEADER_HEIGHT + idx * ROW_HEIGHT;
            return (
              <g key={r.gene}>
                {/* Row label (gene) */}
                <a
                  href={`/gene/${encodeURIComponent(r.gene)}`}
                  onClick={(e) => {
                    if (onRowClick) {
                      e.preventDefault();
                      onRowClick(r.gene);
                    }
                  }}
                >
                  <title>{r.gene} — {r.count}/7 modalities</title>
                  <rect
                    x={0}
                    y={y}
                    width={ROW_LABEL_WIDTH}
                    height={ROW_HEIGHT}
                    fill="transparent"
                  />
                  <text
                    x={8}
                    y={y + ROW_HEIGHT - 5}
                    fontSize={11}
                    fontFamily="monospace"
                    fill="hsl(220 60% 35%)"
                    style={{ textDecoration: "underline" }}
                  >
                    {r.gene}
                  </text>
                </a>
                {/* Modality cells */}
                {r.modalities.map((v, i) => (
                  <rect
                    key={i}
                    x={ROW_LABEL_WIDTH + i * COL_WIDTH + 1}
                    y={y + 1}
                    width={COL_WIDTH - 2}
                    height={ROW_HEIGHT - 2}
                    fill={cellColor(v)}
                    stroke={v > 0.1 ? "hsl(24 60% 55%)" : "hsl(0 0% 92%)"}
                    strokeWidth={0.5}
                  >
                    <title>
                      {r.gene} · M{i + 1}: {v.toFixed(2)}
                    </title>
                  </rect>
                ))}
                {/* Count */}
                <text
                  x={ROW_LABEL_WIDTH + MODALITY_HEADERS.length * COL_WIDTH + 6}
                  y={y + ROW_HEIGHT - 5}
                  fontSize={10}
                  fontFamily="monospace"
                  fill="hsl(0 0% 35%)"
                >
                  {r.count}
                </text>
              </g>
            );
          })}
        </svg>
      </div>
      <p className="text-xs text-muted-foreground">
        Cell intensity: white = absent, dark orange = maximum evidence.
        Click a gene to open its profile; click a column header to filter rows
        where that modality is active.
      </p>
    </div>
  );
}
