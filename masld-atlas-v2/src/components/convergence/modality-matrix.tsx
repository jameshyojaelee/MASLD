"use client";

import { useMemo, useState } from "react";
import { Badge } from "@/components/ui/badge";
import { useChartTooltip, ChartTooltip } from "@/components/charts";
import { sequentialColor } from "@/lib/palette";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
//
// Six MAIN human convergence channels (Fig 5A), in the order emitted by
// scripts/portal/generate_convergence_data.py. Essentiality (DepMap) and
// mouse/cross-species are SUPPLEMENTARY evidence, NOT convergence channels,
// and are intentionally absent here. `tier`/`score`/`rank` are the canonical
// Script-46d convergence values (Tier 1 = genetically validated).

export interface ConvergenceRow {
  gene: string;
  modalities: number[]; // length 6: [bulk, coloc, atac, spatial, singlecell, proteomics]
  count: number; // # of the 6 channels active (> 0.1) — internally consistent with the dots
  tier?: number; // 1 = genetically validated … 4 = weak; 0 = excluded
  tier_label?: string;
  score?: number | null;
  rank?: number | null;
  is_deg: boolean;
  is_coloc: boolean;
  is_druggable: boolean;
}

interface ModalityMatrixProps {
  data: ConvergenceRow[];
  onRowClick?: (gene: string) => void;
}

// Channel headers — order matches `modalities[]`. The six MAIN human modalities.
const MODALITY_HEADERS = [
  { key: "bulk", label: "Bulk", title: "Bulk RNA-seq (disease logFC)" },
  { key: "coloc", label: "COLOC", title: "GWAS–eQTL colocalization (SuSiE PP.H4)" },
  { key: "atac", label: "ATAC", title: "Regulatory / ATAC accessibility" },
  { key: "spatial", label: "Spatial", title: "Spatial (Moran's I / SVG)" },
  { key: "sc", label: "scRNA", title: "Single-cell (cross-cell-type DE)" },
  { key: "prot", label: "Protein", title: "Proteomics (best protein logFC)" },
];
const N_CHANNELS = MODALITY_HEADERS.length; // 6

const ROW_HEIGHT = 18;
const ROW_LABEL_WIDTH = 92;
const COL_WIDTH = 60;
const HEADER_HEIGHT = 48;
const MATRIX_HEIGHT = 580;
const DEFAULT_TOP_N = 500;
const MATRIX_WIDTH = ROW_LABEL_WIDTH + N_CHANNELS * COL_WIDTH + 40;

// Sequential magnitude fill (0 -> light, 1 -> deep). Absent cells use a
// theme-neutral empty color so the ramp reads only as evidence strength.
function cellFill(v: number): string {
  if (v <= 0) return "var(--color-muted)";
  return sequentialColor(Math.min(Math.max(v, 0), 1));
}

type SortMode = "count" | "alpha" | "bulk" | "coloc" | "rank";

type MatrixTip =
  | { kind: "cell"; gene: string; mod: number; v: number }
  | { kind: "row"; gene: string; count: number; tier?: number; tierLabel?: string }
  | { kind: "header"; mod: number };

export function ModalityMatrix({ data, onRowClick }: ModalityMatrixProps) {
  const [query, setQuery] = useState("");
  const [minCount, setMinCount] = useState(0);
  const [requiredModality, setRequiredModality] = useState<number | null>(null);
  const [sortMode, setSortMode] = useState<SortMode>("rank");
  const [showAll, setShowAll] = useState(false);
  const [hover, setHover] = useState<{ gene: string | null; col: number } | null>(
    null
  );
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<MatrixTip>();

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
      case "bulk":
        sorted.sort((a, b) => (b.modalities[0] ?? 0) - (a.modalities[0] ?? 0));
        break;
      case "coloc":
        sorted.sort((a, b) => (b.modalities[1] ?? 0) - (a.modalities[1] ?? 0));
        break;
      case "rank":
        sorted.sort(
          (a, b) =>
            (a.rank ?? Number.MAX_SAFE_INTEGER) - (b.rank ?? Number.MAX_SAFE_INTEGER)
        );
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
    const h = new Array(N_CHANNELS + 1).fill(0);
    data.forEach((r) => {
      h[Math.min(r.count, N_CHANNELS)] += 1;
    });
    return h;
  }, [data]);

  const svgHeight = HEADER_HEIGHT + displayedRows.length * ROW_HEIGHT + 4;

  // Row index of the hovered gene (for the row cross-highlight overlay).
  const hoverRowIdx = useMemo(() => {
    if (!hover?.gene) return -1;
    return displayedRows.findIndex((r) => r.gene === hover.gene);
  }, [hover, displayedRows]);

  function clearHover() {
    setHover(null);
    hide();
  }

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
          genes across{" "}
          <span className="font-mono font-semibold">{N_CHANNELS}</span> main human
          channels. Highest convergence:{" "}
          <span className="font-mono font-semibold">{maxCount}</span>/{N_CHANNELS}.
        </p>
        <p className="mt-1 text-xs text-muted-foreground">
          By channel count:{" "}
          {countHistogram
            .map((n, i) => (n > 0 ? `${i}/${N_CHANNELS}=${n.toLocaleString()}` : null))
            .filter(Boolean)
            .join(" · ")}
          . Essentiality and cross-species are supplementary, not convergence
          channels.
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
          <span className="text-xs text-muted-foreground">Min channels:</span>
          <input
            type="range"
            min={0}
            max={N_CHANNELS}
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
            <option value="rank">Convergence rank</option>
            <option value="count">Channel count</option>
            <option value="alpha">Gene (A–Z)</option>
            <option value="bulk">Bulk DE</option>
            <option value="coloc">Genetic (COLOC)</option>
          </select>
        </label>
        {requiredModality !== null && (
          <Badge variant="outline" className="text-[10px]">
            {MODALITY_HEADERS[requiredModality].label} required
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

      {/* Matrix (tooltip anchored to this non-scrolling wrapper). */}
      <div ref={wrapperRef} className="relative">
        <div
          className="relative overflow-auto rounded-md border border-border bg-background"
          style={{ maxHeight: MATRIX_HEIGHT }}
          onMouseLeave={clearHover}
        >
          <svg width={MATRIX_WIDTH} height={svgHeight} style={{ display: "block" }}>
            {/* Column cross-highlight */}
            {hover && hover.col >= 0 && (
              <rect
                x={ROW_LABEL_WIDTH + hover.col * COL_WIDTH}
                y={HEADER_HEIGHT}
                width={COL_WIDTH}
                height={displayedRows.length * ROW_HEIGHT}
                fill="var(--color-accent)"
                opacity={0.5}
              />
            )}
            {/* Row cross-highlight */}
            {hoverRowIdx >= 0 && (
              <rect
                x={0}
                y={HEADER_HEIGHT + hoverRowIdx * ROW_HEIGHT}
                width={MATRIX_WIDTH}
                height={ROW_HEIGHT}
                fill="var(--color-accent)"
                opacity={0.5}
              />
            )}

            {/* Header row */}
            <g>
              <rect
                x={0}
                y={0}
                width={MATRIX_WIDTH}
                height={HEADER_HEIGHT}
                fill="var(--color-card)"
                stroke="var(--color-border)"
              />
              {MODALITY_HEADERS.map((m, i) => {
                const active = requiredModality === i;
                return (
                  <g
                    key={m.key}
                    transform={`translate(${ROW_LABEL_WIDTH + i * COL_WIDTH},0)`}
                    onClick={() => setRequiredModality(active ? null : i)}
                    onMouseMove={(e) => {
                      setHover({ gene: null, col: i });
                      show(e, { kind: "header", mod: i });
                    }}
                    style={{ cursor: "pointer" }}
                  >
                    <rect
                      x={1}
                      y={1}
                      width={COL_WIDTH - 2}
                      height={HEADER_HEIGHT - 2}
                      fill={active ? "var(--color-accent)" : "transparent"}
                    />
                    <text
                      x={COL_WIDTH / 2}
                      y={HEADER_HEIGHT / 2 + 4}
                      textAnchor="middle"
                      fontSize={10}
                      fontFamily="var(--font-mono)"
                      fontWeight={600}
                      fill="var(--color-foreground)"
                    >
                      {m.label}
                    </text>
                  </g>
                );
              })}
              <text
                x={ROW_LABEL_WIDTH + N_CHANNELS * COL_WIDTH + 6}
                y={HEADER_HEIGHT / 2 + 4}
                fontSize={10}
                fill="var(--color-muted-foreground)"
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
                    href={`#/gene?symbol=${encodeURIComponent(r.gene)}`}
                    onClick={(e) => {
                      if (onRowClick) {
                        e.preventDefault();
                        onRowClick(r.gene);
                      }
                    }}
                    onMouseMove={(e) => {
                      setHover({ gene: r.gene, col: -1 });
                      show(e, {
                        kind: "row",
                        gene: r.gene,
                        count: r.count,
                        tier: r.tier,
                        tierLabel: r.tier_label,
                      });
                    }}
                  >
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
                      fontFamily="var(--font-mono)"
                      fill="var(--color-primary)"
                      style={{ textDecoration: "underline" }}
                    >
                      {r.gene}
                    </text>
                  </a>
                  {/* Channel cells */}
                  {r.modalities.slice(0, N_CHANNELS).map((v, i) => (
                    <rect
                      key={i}
                      x={ROW_LABEL_WIDTH + i * COL_WIDTH + 1}
                      y={y + 1}
                      width={COL_WIDTH - 2}
                      height={ROW_HEIGHT - 2}
                      fill={cellFill(v)}
                      stroke="var(--color-border)"
                      strokeWidth={0.5}
                      strokeOpacity={v > 0.1 ? 0.9 : 0.4}
                      onMouseMove={(e) => {
                        setHover({ gene: r.gene, col: i });
                        show(e, { kind: "cell", gene: r.gene, mod: i, v });
                      }}
                    />
                  ))}
                  {/* Count */}
                  <text
                    x={ROW_LABEL_WIDTH + N_CHANNELS * COL_WIDTH + 6}
                    y={y + ROW_HEIGHT - 5}
                    fontSize={10}
                    fontFamily="var(--font-mono)"
                    fill="var(--color-muted-foreground)"
                  >
                    {r.count}
                  </text>
                </g>
              );
            })}
          </svg>
        </div>

        {tooltip && (
          <ChartTooltip left={tooltip.left} top={tooltip.top}>
            {tooltip.data.kind === "cell" && (
              <>
                <div className="font-medium">
                  {tooltip.data.gene} · {MODALITY_HEADERS[tooltip.data.mod].label}
                </div>
                <div className="text-muted-foreground">
                  {MODALITY_HEADERS[tooltip.data.mod].title}
                </div>
                <div>evidence {tooltip.data.v.toFixed(2)}</div>
              </>
            )}
            {tooltip.data.kind === "row" && (
              <>
                <div className="font-medium">
                  {tooltip.data.gene} — {tooltip.data.count}/{N_CHANNELS} channels
                </div>
                {tooltip.data.tierLabel && (
                  <div className="text-muted-foreground">{tooltip.data.tierLabel}</div>
                )}
              </>
            )}
            {tooltip.data.kind === "header" && (
              <>
                <div className="font-medium">
                  {MODALITY_HEADERS[tooltip.data.mod].label}
                </div>
                <div className="text-muted-foreground">
                  {MODALITY_HEADERS[tooltip.data.mod].title}
                </div>
                <div className="text-muted-foreground">
                  Click to filter rows where this channel is active.
                </div>
              </>
            )}
          </ChartTooltip>
        )}
      </div>
      <p className="text-xs text-muted-foreground">
        Six main human channels: bulk RNA-seq, GWAS–eQTL colocalization,
        regulatory/ATAC, spatial, single-cell, proteomics. Cell intensity: light =
        absent, deep = maximum evidence. Hover a cell to cross-highlight its gene
        row and channel column; click a gene to open its profile; click a column
        header to filter.
      </p>
    </div>
  );
}
