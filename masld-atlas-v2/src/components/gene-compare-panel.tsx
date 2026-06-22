"use client";

import { useEffect, useState, useMemo } from "react";
import { useAppStore } from "@/lib/store";
import { getGeneBySymbol } from "@/lib/search-index";
import { EVIDENCE_SOURCES } from "@/lib/colors";
import type { GeneIndexEntry, EvidenceStrengths } from "@/lib/types";

/** Up to 5 distinct colors for compared genes */
const COMPARE_COLORS = ["#3b82f6", "#ef4444", "#22c55e", "#f59e0b", "#8b5cf6"];

/**
 * gene_index.json uses the M1-M7 aligned key scheme (s2_mouse, s3_genetic,
 * s4_essential, s5_epigenomic, s6_spatial, s7_singlecell). The TypeScript
 * EvidenceStrengths keeps the legacy keys (s2_genetic, s3_essential, ...).
 * Map back so radar/table lookups find the right values.
 */
function normalizeEvidence(e: Record<string, number> | undefined): EvidenceStrengths {
  const raw = e ?? {};
  return {
    s1_human: raw.s1_human ?? 0,
    s2_genetic: raw.s3_genetic ?? raw.s2_genetic ?? 0,
    s3_essential: raw.s4_essential ?? raw.s3_essential ?? 0,
    s4_epigenomic: raw.s5_epigenomic ?? raw.s4_epigenomic ?? 0,
    s5_spatial: raw.s6_spatial ?? raw.s5_spatial ?? 0,
    s6_singlecell: raw.s7_singlecell ?? raw.s6_singlecell ?? 0,
    s7_mouse: raw.s2_mouse ?? raw.s7_mouse ?? 0,
    s8_proteomics: raw.s8_proteomics ?? 0,
  };
}

export function GeneComparePanel() {
  const { compareGenes, removeCompareGene, clearCompareGenes } = useAppStore();
  const [geneData, setGeneData] = useState<GeneIndexEntry[]>([]);
  const [collapsed, setCollapsed] = useState(false);

  useEffect(() => {
    if (compareGenes.length < 2) {
      setGeneData([]);
      return;
    }
    let cancelled = false;
    Promise.all(compareGenes.map((s) => getGeneBySymbol(s))).then((results) => {
      if (!cancelled) {
        setGeneData(results.filter(Boolean) as GeneIndexEntry[]);
      }
    });
    return () => {
      cancelled = true;
    };
  }, [compareGenes]);

  if (compareGenes.length < 2) return null;

  return (
    <div
      className="fixed bottom-0 left-0 right-0 z-40 border-t bg-background/95 backdrop-blur-sm shadow-lg transition-transform duration-200"
      style={{ transform: collapsed ? "translateY(calc(100% - 40px))" : "translateY(0)" }}
    >
      {/* Header bar */}
      <div className="flex items-center justify-between border-b px-4 py-2">
        <div className="flex items-center gap-2">
          <button
            onClick={() => setCollapsed(!collapsed)}
            className="inline-flex items-center gap-1.5 text-sm font-medium hover:text-foreground text-muted-foreground transition-colors"
            aria-label={collapsed ? "Expand comparison panel" : "Collapse comparison panel"}
          >
            <svg
              className="h-4 w-4 transition-transform"
              style={{ transform: collapsed ? "rotate(180deg)" : "rotate(0deg)" }}
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
              strokeWidth={2}
            >
              <path strokeLinecap="round" strokeLinejoin="round" d="M19 9l-7 7-7-7" />
            </svg>
            Compare ({compareGenes.length})
          </button>
        </div>

        {/* Gene tabs */}
        <div className="flex items-center gap-1.5">
          {compareGenes.map((symbol, i) => (
            <span
              key={symbol}
              className="inline-flex items-center gap-1 rounded-full border px-2.5 py-0.5 text-xs font-medium"
              style={{ borderColor: COMPARE_COLORS[i % COMPARE_COLORS.length] }}
            >
              <span
                className="inline-block h-2 w-2 rounded-full"
                style={{ backgroundColor: COMPARE_COLORS[i % COMPARE_COLORS.length] }}
              />
              {symbol}
              <button
                onClick={() => removeCompareGene(symbol)}
                className="ml-0.5 text-muted-foreground hover:text-foreground"
                aria-label={`Remove ${symbol}`}
              >
                &times;
              </button>
            </span>
          ))}
          <button
            onClick={clearCompareGenes}
            className="ml-2 text-xs text-muted-foreground hover:text-foreground transition-colors"
          >
            Clear all
          </button>
        </div>
      </div>

      {/* Main content */}
      {!collapsed && (
        <div className="flex gap-6 overflow-x-auto px-4 py-4" style={{ maxHeight: 280 }}>
          {/* Overlaid radar chart */}
          <OverlaidRadar geneData={geneData} compareGenes={compareGenes} />

          {/* Comparison table */}
          <ComparisonTable geneData={geneData} compareGenes={compareGenes} />
        </div>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Overlaid Radar
// ---------------------------------------------------------------------------

function OverlaidRadar({
  geneData,
  compareGenes,
}: {
  geneData: GeneIndexEntry[];
  compareGenes: string[];
}) {
  const size = 220;
  const cx = size / 2;
  const cy = size / 2;
  const radius = size / 2 - 24;
  const axes = EVIDENCE_SOURCES;
  const n = axes.length;

  const gridLevels = [0.25, 0.5, 0.75, 1];

  /** Compute polygon points string for one gene's evidence values */
  const computePolygon = useMemo(() => {
    return (evidence: EvidenceStrengths) =>
      axes
        .map((source, i) => {
          const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
          const val = evidence[source.key] ?? 0;
          const r = radius * val;
          return `${cx + r * Math.cos(angle)},${cy + r * Math.sin(angle)}`;
        })
        .join(" ");
  }, [axes, n, cx, cy, radius]);

  /** Axis label positions (slightly beyond radius) */
  const labels = useMemo(
    () =>
      axes.map((source, i) => {
        const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
        const labelR = radius + 16;
        return {
          x: cx + labelR * Math.cos(angle),
          y: cy + labelR * Math.sin(angle),
          label: source.short,
        };
      }),
    [axes, n, cx, cy, radius],
  );

  return (
    <div className="flex-shrink-0">
      <svg width={size} height={size} viewBox={`0 0 ${size} ${size}`}>
        {/* Grid rings */}
        {gridLevels.map((level) => (
          <polygon
            key={level}
            points={axes
              .map((_, i) => {
                const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
                const r = radius * level;
                return `${cx + r * Math.cos(angle)},${cy + r * Math.sin(angle)}`;
              })
              .join(" ")}
            fill="none"
            stroke="currentColor"
            strokeWidth={0.5}
            opacity={0.1}
          />
        ))}

        {/* Axis lines */}
        {axes.map((_, i) => {
          const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
          return (
            <line
              key={i}
              x1={cx}
              y1={cy}
              x2={cx + radius * Math.cos(angle)}
              y2={cy + radius * Math.sin(angle)}
              stroke="currentColor"
              strokeWidth={0.5}
              opacity={0.1}
            />
          );
        })}

        {/* One polygon per gene */}
        {geneData.map((gene, idx) => {
          const colorIdx = compareGenes.indexOf(gene.symbol);
          const color = COMPARE_COLORS[(colorIdx >= 0 ? colorIdx : idx) % COMPARE_COLORS.length];
          return (
            <polygon
              key={gene.symbol}
              points={computePolygon(normalizeEvidence(gene.evidence as unknown as Record<string, number>))}
              fill={color}
              fillOpacity={0.15}
              stroke={color}
              strokeWidth={1.5}
              strokeOpacity={0.7}
            />
          );
        })}

        {/* Axis labels */}
        {labels.map((l) => (
          <text
            key={l.label}
            x={l.x}
            y={l.y}
            textAnchor="middle"
            dominantBaseline="central"
            className="fill-muted-foreground"
            fontSize={9}
          >
            {l.label}
          </text>
        ))}
      </svg>

      {/* Legend */}
      <div className="mt-1 flex flex-wrap gap-x-3 gap-y-0.5 px-2">
        {geneData.map((gene, idx) => {
          const colorIdx = compareGenes.indexOf(gene.symbol);
          const color = COMPARE_COLORS[(colorIdx >= 0 ? colorIdx : idx) % COMPARE_COLORS.length];
          return (
            <span key={gene.symbol} className="inline-flex items-center gap-1 text-[10px]">
              <span
                className="inline-block h-2 w-2 rounded-full"
                style={{ backgroundColor: color }}
              />
              {gene.symbol}
            </span>
          );
        })}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Comparison Table
// ---------------------------------------------------------------------------

function ComparisonTable({
  geneData,
  compareGenes,
}: {
  geneData: GeneIndexEntry[];
  compareGenes: string[];
}) {
  /** Ordered gene data matching compareGenes order */
  const ordered = useMemo(() => {
    return compareGenes
      .map((s) => geneData.find((g) => g.symbol === s))
      .filter(Boolean) as GeneIndexEntry[];
  }, [compareGenes, geneData]);

  if (ordered.length === 0) return null;

  const rows: { label: string; render: (gene: GeneIndexEntry) => string }[] = [
    {
      label: "logFC",
      render: (g) => (g.bulk_logfc != null ? g.bulk_logfc.toFixed(2) : "N/A"),
    },
    {
      label: "padj",
      render: (g) =>
        g.bulk_padj != null
          ? g.bulk_padj < 0.001
            ? g.bulk_padj.toExponential(1)
            : g.bulk_padj.toFixed(3)
          : "N/A",
    },
    {
      label: "Modalities",
      render: (g) => `${g.layers_active}/7`,
    },
    {
      label: "DEG",
      render: (g) => (g.is_deg ? "Yes" : "No"),
    },
    {
      label: "Conserved Core",
      render: (g) => (g.is_conserved_core ? "Yes" : "No"),
    },
    {
      label: "Druggable",
      render: (g) => (g.dgidb_druggable ? "Yes" : "No"),
    },
    {
      label: "Biotype",
      render: (g) => g.biotype ?? "unknown",
    },
    {
      label: "Sex Class",
      render: (g) => g.sex_class ?? "N/A",
    },
  ];

  return (
    <div className="flex-1 overflow-x-auto">
      <table className="w-full text-xs">
        <thead>
          <tr className="border-b">
            <th className="py-1.5 pr-4 text-left font-medium text-muted-foreground">Metric</th>
            {ordered.map((gene, idx) => {
              const color = COMPARE_COLORS[idx % COMPARE_COLORS.length];
              return (
                <th
                  key={gene.symbol}
                  className="py-1.5 px-3 text-left font-medium"
                  style={{ color }}
                >
                  {gene.symbol}
                </th>
              );
            })}
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => (
            <tr key={row.label} className="border-b border-border/50">
              <td className="py-1.5 pr-4 text-muted-foreground">{row.label}</td>
              {ordered.map((gene) => (
                <td key={gene.symbol} className="py-1.5 px-3 font-mono">
                  {row.render(gene)}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
