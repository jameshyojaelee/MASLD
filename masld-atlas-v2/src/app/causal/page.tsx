"use client";

// LEGACY: this route is absorbed by /genetics — retire after migration. Kept
// live for now; changes here are correctness-only (no expanded surface).

import { useEffect, useState, useMemo } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { PageContainer } from "@/components/page-container";
import { PageHeader } from "@/components/page-header";
import { SkeletonBlock } from "@/components/states";
import { useChartTooltip, ChartTooltip } from "@/components/charts";
import { dataUrl } from "@/lib/data-base";
import { COLOC_SUSIE, GWAS_COUNT, fmt } from "@/lib/atlas-constants";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface MotifDisrupted {
  tf: string;
  effect: string;
  in_regulon: boolean;
  alleleDiff: number;
}

interface Variant {
  id: string;
  chr: string;
  pos: number;
  pip: number;
  rec_pip: number;
  cell_types: string[];
  nearest_gene: string;
  distance: number;
  motifs_disrupted: MotifDisrupted[];
  coloc_pp4: number;
  coloc_gwas: string;
}

interface CellTypeEnrichment {
  cell_type: string;
  count: number;
  enrichment: number;
  padj: number;
}

interface GwasAtacData {
  variants: Variant[];
  summary: {
    total_variants: number;
    in_peaks: number;
    with_motif_disruption: number;
    disease_regulon_tfs: number;
  };
  cell_type_enrichment: CellTypeEnrichment[];
}

// ---------------------------------------------------------------------------
// Hardcoded framework stats
// ---------------------------------------------------------------------------

interface FrameworkStat {
  value: string;
  label: string;
  sublabel: string;
}

const FRAMEWORK_STATS: FrameworkStat[] = [
  {
    value: fmt(COLOC_SUSIE),
    label: "SuSiE effector genes",
    sublabel: "Tier-1/2 liver-specific, PP.H4 > 0.5",
  },
  {
    value: String(GWAS_COUNT),
    label: "GWAS studies",
    sublabel: "EUR, AFR, EAS, AMR, SAS ancestries",
  },
  {
    value: "5",
    label: "Ancestries",
    sublabel: "European, African, East Asian, Admixed American, South Asian",
  },
  {
    value: "248",
    label: "Cross-ancestry genes",
    sublabel: "PP.H4 replicated in ≥2 ancestries",
  },
];

// ---------------------------------------------------------------------------
// Top regulatory TF findings
// ---------------------------------------------------------------------------

interface TopTf {
  tf: string;
  variants: number;
  description: string;
}

const TOP_TFS: TopTf[] = [
  {
    tf: "HNF4A",
    variants: 12,
    description:
      "Master hepatocyte TF; motif disrupted by 12 GWAS variants in scATAC peaks",
  },
  {
    tf: "RORA",
    variants: 11,
    description:
      "Circadian/metabolic regulator; 11 variants disrupt binding across multiple cell types",
  },
  {
    tf: "THRB",
    variants: 8,
    description:
      "Thyroid hormone receptor (resmetirom target); 8 variants in regulatory elements",
  },
];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function formatPadj(p: number): string {
  if (p < 0.001) return p.toExponential(1);
  return p.toFixed(3);
}

function formatCellType(ct: string): string {
  return ct.replace(/_/g, " ");
}

const ROWS_PER_PAGE = 20;

// ---------------------------------------------------------------------------
// Cell Type Enrichment Bar Chart (SVG)
// ---------------------------------------------------------------------------

const BAR_CHART_W = 700;
const BAR_CHART_H_PER_ROW = 36;
const BAR_MARGIN = { top: 10, right: 180, bottom: 10, left: 130 };

function CellTypeEnrichmentChart({
  data,
}: {
  data: CellTypeEnrichment[];
}) {
  const { wrapperRef, tooltip, show, hide } =
    useChartTooltip<CellTypeEnrichment>();

  // Sort by count descending
  const sorted = [...data].sort((a, b) => b.count - a.count);
  const maxCount = Math.max(...sorted.map((d) => d.count));

  const innerW = BAR_CHART_W - BAR_MARGIN.left - BAR_MARGIN.right;
  const totalH =
    BAR_MARGIN.top + sorted.length * BAR_CHART_H_PER_ROW + BAR_MARGIN.bottom;

  return (
    <div ref={wrapperRef} className="relative">
      <svg
        width={BAR_CHART_W}
        height={totalH}
        className="max-w-full overflow-visible"
        style={{ fontFamily: "inherit" }}
      >
        <g transform={`translate(${BAR_MARGIN.left},${BAR_MARGIN.top})`}>
          {sorted.map((d, i) => {
            const barW = (d.count / maxCount) * innerW;
            const y = i * BAR_CHART_H_PER_ROW;
            const barH = BAR_CHART_H_PER_ROW - 6;
            const isSig = d.padj < 0.05;

            return (
              <g
                key={d.cell_type}
                onMouseMove={(e) => show(e, d)}
                onMouseLeave={hide}
              >
                {/* Label */}
                <text
                  x={-8}
                  y={y + barH / 2 + 4}
                  textAnchor="end"
                  fontSize={11}
                  fill="var(--color-foreground)"
                  fontWeight={500}
                >
                  {formatCellType(d.cell_type)}
                </text>
                {/* Bar: significance by palette hue (genetic modality) + opacity;
                    non-significant falls back to control gray. */}
                <rect
                  x={0}
                  y={y}
                  width={barW}
                  height={barH}
                  rx={3}
                  fill={isSig ? "var(--color-s2-genetic)" : "var(--color-control)"}
                  opacity={isSig ? 0.85 : 0.4}
                />
                {/* Count inside bar */}
                <text
                  x={barW + 6}
                  y={y + barH / 2 + 4}
                  fontSize={10}
                  fill="var(--color-foreground)"
                  fontWeight={600}
                >
                  {d.count}
                </text>
                {/* OR + padj annotation */}
                <text
                  x={barW + 40}
                  y={y + barH / 2 + 4}
                  fontSize={10}
                  fill="var(--color-muted-foreground)"
                >
                  OR={d.enrichment.toFixed(2)}{" "}
                  {isSig ? `padj=${formatPadj(d.padj)}` : "ns"}
                </text>
              </g>
            );
          })}
        </g>
      </svg>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium">
            {formatCellType(tooltip.data.cell_type)}
          </div>
          <div>{tooltip.data.count} variants</div>
          <div>
            OR {tooltip.data.enrichment.toFixed(2)} ·{" "}
            {tooltip.data.padj < 0.05
              ? `padj ${formatPadj(tooltip.data.padj)}`
              : "n.s."}
          </div>
        </ChartTooltip>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// PIP Bar (inline)
// ---------------------------------------------------------------------------

function PipBar({ pip }: { pip: number }) {
  const widthPct = Math.max(pip * 100, 2);
  return (
    <div className="flex items-center gap-2">
      <div className="h-2 w-16 overflow-hidden rounded-full bg-muted">
        <div
          className="h-full rounded-full bg-primary"
          style={{ width: `${widthPct}%` }}
        />
      </div>
      <span className="font-mono text-xs">{pip.toFixed(3)}</span>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Motif Detail (expanded row)
// ---------------------------------------------------------------------------

function MotifDetail({ motifs }: { motifs: MotifDisrupted[] }) {
  if (motifs.length === 0) {
    return (
      <span className="text-xs text-muted-foreground italic">
        No motifs disrupted
      </span>
    );
  }

  // Deduplicate by TF name (data may have duplicates)
  const seen = new Set<string>();
  const unique: MotifDisrupted[] = [];
  for (const m of motifs) {
    if (!seen.has(m.tf)) {
      seen.add(m.tf);
      unique.push(m);
    }
  }

  return (
    <div className="space-y-2">
      {unique.map((m, i) => (
        <div
          key={`${m.tf}-${i}`}
          className="flex items-center gap-3 text-xs"
        >
          <Badge
            variant={m.in_regulon ? "default" : "outline"}
            className="text-[10px]"
          >
            {m.tf}
          </Badge>
          <span className="text-muted-foreground">
            Effect:{" "}
            <span
              className={
                m.effect === "strong" ? "font-semibold text-foreground" : ""
              }
            >
              {m.effect}
            </span>
          </span>
          <span className="text-muted-foreground">
            alleleDiff: {m.alleleDiff.toFixed(3)}
          </span>
          {m.in_regulon && (
            <Badge variant="destructive" className="text-[10px]">
              Disease regulon
            </Badge>
          )}
        </div>
      ))}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function CausalPage() {
  const [data, setData] = useState<GwasAtacData | null>(null);
  const [loading, setLoading] = useState(true);
  const [cellTypeFilter, setCellTypeFilter] = useState<string>("all");
  const [page, setPage] = useState(0);
  const [expandedRow, setExpandedRow] = useState<string | null>(null);

  useEffect(() => {
    fetch(dataUrl("gwas_atac_browser.json"))
      .then((res) => res.json())
      .then((json: GwasAtacData) => {
        setData(json);
        setLoading(false);
      });
  }, []);

  // All unique cell types from variants
  const cellTypes = useMemo(() => {
    if (!data) return [];
    const s = new Set<string>();
    for (const v of data.variants) {
      for (const ct of v.cell_types) s.add(ct);
    }
    return Array.from(s).sort();
  }, [data]);

  // Filtered + sorted variants
  const filteredVariants = useMemo(() => {
    if (!data) return [];
    let filtered = data.variants;
    if (cellTypeFilter !== "all") {
      filtered = filtered.filter((v) =>
        v.cell_types.includes(cellTypeFilter)
      );
    }
    return [...filtered].sort((a, b) => b.pip - a.pip);
  }, [data, cellTypeFilter]);

  // Pagination
  const totalPages = Math.ceil(filteredVariants.length / ROWS_PER_PAGE);
  const pagedVariants = filteredVariants.slice(
    page * ROWS_PER_PAGE,
    (page + 1) * ROWS_PER_PAGE
  );

  // Reset page when filter changes
  useEffect(() => {
    setPage(0);
    setExpandedRow(null);
  }, [cellTypeFilter]);

  return (
    <PageContainer>
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <PageHeader
        title="Causal Architecture"
        description={
          <>
            COLOC + TWAS across {GWAS_COUNT} GWAS (5 ancestries), multi-ancestry
            replication, regulatory variant mapping
          </>
        }
      />

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Causal Framework Overview                                 */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Causal Framework Overview
        </h2>
        <p className="mb-6 text-sm leading-relaxed text-muted-foreground">
          We integrate colocalization (COLOC), transcriptome-wide association
          (TWAS), and fine-mapping to identify genetically causal MASLD genes.
          SuSiE fine-mapping with ABF COLOC fallback was applied across{" "}
          {GWAS_COUNT} GWAS spanning 5 ancestries (European, African, East
          Asian, Admixed American, South Asian) for cross-ancestry
          replication. Regulatory variant mapping intersects credible set
          variants with single-cell ATAC-seq peaks to identify causal
          non-coding mechanisms.
        </p>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
          {FRAMEWORK_STATS.map((stat) => (
            <div
              key={stat.label}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm"
            >
              <p className="font-mono text-2xl font-bold text-primary">
                {stat.value}
              </p>
              <p className="mt-0.5 text-sm font-semibold">{stat.label}</p>
              <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                {stat.sublabel}
              </p>
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 2: GWAS-ATAC Regulatory Variant Browser                      */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          GWAS-ATAC Regulatory Variant Browser
        </h2>

        {loading ? (
          <div className="space-y-6">
            <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
              {Array.from({ length: 4 }).map((_, i) => (
                <SkeletonBlock key={i} className="h-24 rounded-lg" />
              ))}
            </div>
            <SkeletonBlock className="h-56 w-full rounded-lg" />
            <SkeletonBlock className="h-72 w-full rounded-lg" />
          </div>
        ) : data ? (
          <>
            {/* Summary cards */}
            <div className="mb-6 grid grid-cols-2 gap-4 sm:grid-cols-4">
              <div className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm">
                <p className="font-mono text-2xl font-bold text-primary">
                  {data.summary.total_variants.toLocaleString()}
                </p>
                <p className="mt-0.5 text-sm font-semibold">
                  Credible set variants
                </p>
                <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                  SuSiE + CARMA 95% CS
                </p>
              </div>
              <div className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm">
                <p className="font-mono text-2xl font-bold text-primary">
                  {data.summary.in_peaks.toLocaleString()}
                </p>
                <p className="mt-0.5 text-sm font-semibold">
                  In scATAC peaks
                </p>
                <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                  {((data.summary.in_peaks / data.summary.total_variants) * 100).toFixed(1)}% of credible set
                </p>
              </div>
              <div className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm">
                <p className="font-mono text-2xl font-bold text-primary">
                  {data.summary.with_motif_disruption}
                </p>
                <p className="mt-0.5 text-sm font-semibold">
                  Motif disruption
                </p>
                <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                  motifbreakR + FIMO validated
                </p>
              </div>
              <div className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm">
                <p className="font-mono text-2xl font-bold text-primary">
                  {data.summary.disease_regulon_tfs}
                </p>
                <p className="mt-0.5 text-sm font-semibold">
                  Disease regulon TFs
                </p>
                <p className="mt-1 text-[11px] leading-tight text-muted-foreground">
                  SCENIC+ regulon overlap
                </p>
              </div>
            </div>

            {/* Cell type enrichment chart */}
            <div className="mb-6">
              <h3 className="mb-3 text-sm font-semibold text-muted-foreground uppercase tracking-wider">
                Cell-Type Enrichment
              </h3>
              <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
                <CellTypeEnrichmentChart data={data.cell_type_enrichment} />
                <p className="mt-3 text-xs text-muted-foreground">
                  Horizontal bars show variant count per cell type. Blue bars are
                  significantly enriched (padj &lt; 0.05). OR = enrichment odds
                  ratio vs genomic background.
                </p>
              </div>
            </div>

            {/* Variant table */}
            <div className="mb-4">
              <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
                <h3 className="text-sm font-semibold text-muted-foreground uppercase tracking-wider">
                  Variants in scATAC Peaks
                </h3>
                <div className="flex items-center gap-2">
                  <label
                    htmlFor="ct-filter"
                    className="text-xs text-muted-foreground"
                  >
                    Filter by cell type:
                  </label>
                  <select
                    id="ct-filter"
                    value={cellTypeFilter}
                    onChange={(e) => setCellTypeFilter(e.target.value)}
                    className="rounded-md border border-border bg-card px-2 py-1 text-xs focus:outline-none focus:ring-1 focus:ring-ring"
                  >
                    <option value="all">All cell types</option>
                    {cellTypes.map((ct) => (
                      <option key={ct} value={ct}>
                        {formatCellType(ct)}
                      </option>
                    ))}
                  </select>
                </div>
              </div>

              <div className="overflow-x-auto rounded-lg border border-border">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-border bg-muted/50">
                      {[
                        "Variant ID",
                        "Chr",
                        "Position",
                        "PIP",
                        "Cell Type(s)",
                        "Nearest Gene",
                        "Dist",
                        "Motifs",
                      ].map((h) => (
                        <th
                          key={h}
                          className="px-3 py-2.5 text-left text-xs font-medium text-muted-foreground"
                        >
                          {h}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {pagedVariants.map((v) => {
                      const isExpanded = expandedRow === v.id;
                      // Deduplicate motif TFs for badge display
                      const uniqueTfs = new Map<string, MotifDisrupted>();
                      for (const m of v.motifs_disrupted) {
                        if (!uniqueTfs.has(m.tf)) uniqueTfs.set(m.tf, m);
                      }
                      const tfList = Array.from(uniqueTfs.values());

                      return (
                        <tr key={v.id} className="group">
                          <td colSpan={8} className="p-0">
                            {/* Main row */}
                            <button
                              type="button"
                              onClick={() =>
                                setExpandedRow(isExpanded ? null : v.id)
                              }
                              className="flex w-full items-center border-b border-border/50 transition-colors hover:bg-muted/30"
                            >
                              <span className="w-[160px] shrink-0 px-3 py-2.5 text-left font-mono text-xs">
                                {v.id}
                              </span>
                              <span className="w-[56px] shrink-0 px-3 py-2.5 text-left font-mono text-xs">
                                {v.chr.replace("chr", "")}
                              </span>
                              <span className="w-[100px] shrink-0 px-3 py-2.5 text-left font-mono text-xs">
                                {v.pos.toLocaleString()}
                              </span>
                              <span className="w-[130px] shrink-0 px-3 py-2.5">
                                <PipBar pip={v.pip} />
                              </span>
                              <span className="w-[140px] shrink-0 px-3 py-2.5 text-left text-xs">
                                <span className="flex flex-wrap gap-1">
                                  {v.cell_types.length <= 2 ? (
                                    v.cell_types.map((ct) => (
                                      <Badge
                                        key={ct}
                                        variant="outline"
                                        className="text-[9px]"
                                      >
                                        {formatCellType(ct)}
                                      </Badge>
                                    ))
                                  ) : (
                                    <Badge
                                      variant="outline"
                                      className="text-[9px]"
                                    >
                                      {v.cell_types.length} cell types
                                    </Badge>
                                  )}
                                </span>
                              </span>
                              <span className="w-[100px] shrink-0 px-3 py-2.5 text-left text-xs font-semibold text-primary">
                                <Link
                                  href={`/gene?symbol=${encodeURIComponent(v.nearest_gene)}`}
                                  onClick={(e) => e.stopPropagation()}
                                  className="hover:underline"
                                >
                                  {v.nearest_gene}
                                </Link>
                              </span>
                              <span className="w-[60px] shrink-0 px-3 py-2.5 text-right font-mono text-xs text-muted-foreground">
                                {v.distance === 0
                                  ? "0"
                                  : v.distance.toLocaleString()}
                              </span>
                              <span className="flex-1 px-3 py-2.5 text-left">
                                <span className="flex flex-wrap gap-1">
                                  {tfList.length === 0 ? (
                                    <span className="text-xs text-muted-foreground">
                                      --
                                    </span>
                                  ) : (
                                    tfList.map((m) => (
                                      <Badge
                                        key={m.tf}
                                        variant={
                                          m.in_regulon
                                            ? "destructive"
                                            : "secondary"
                                        }
                                        className="text-[9px]"
                                      >
                                        {m.tf}
                                      </Badge>
                                    ))
                                  )}
                                </span>
                              </span>
                            </button>
                            {/* Expanded detail */}
                            {isExpanded && (
                              <div className="border-b border-border bg-muted/20 px-4 py-3">
                                <div className="mb-2 flex flex-wrap gap-4 text-xs text-muted-foreground">
                                  <span>
                                    COLOC PP.H4:{" "}
                                    <span className="font-mono font-semibold text-foreground">
                                      {v.coloc_pp4.toFixed(4)}
                                    </span>
                                  </span>
                                  <span>
                                    GWAS:{" "}
                                    <span className="font-mono text-foreground">
                                      {v.coloc_gwas}
                                    </span>
                                  </span>
                                  <span>
                                    Recommended PIP:{" "}
                                    <span className="font-mono font-semibold text-foreground">
                                      {v.rec_pip.toFixed(3)}
                                    </span>
                                  </span>
                                  <span>
                                    Cell types:{" "}
                                    <span className="text-foreground">
                                      {v.cell_types
                                        .map(formatCellType)
                                        .join(", ")}
                                    </span>
                                  </span>
                                </div>
                                <div className="text-xs font-semibold text-muted-foreground mb-1.5">
                                  Motif Disruptions
                                </div>
                                <MotifDetail motifs={v.motifs_disrupted} />
                              </div>
                            )}
                          </td>
                        </tr>
                      );
                    })}
                    {pagedVariants.length === 0 && (
                      <tr>
                        <td
                          colSpan={8}
                          className="px-4 py-8 text-center text-sm text-muted-foreground"
                        >
                          No variants found for this filter.
                        </td>
                      </tr>
                    )}
                  </tbody>
                </table>
              </div>

              {/* Pagination */}
              <div className="mt-3 flex items-center justify-between">
                <span className="text-xs text-muted-foreground">
                  {filteredVariants.length} variant{filteredVariants.length !== 1 ? "s" : ""}{" "}
                  {cellTypeFilter !== "all" &&
                    `in ${formatCellType(cellTypeFilter)}`}
                </span>
                <div className="flex items-center gap-2">
                  <button
                    type="button"
                    disabled={page === 0}
                    onClick={() => setPage((p) => p - 1)}
                    className="rounded-md border border-border px-3 py-1 text-xs font-medium transition-colors hover:bg-muted disabled:opacity-40 disabled:cursor-not-allowed"
                  >
                    Previous
                  </button>
                  <span className="text-xs text-muted-foreground">
                    Page {page + 1} of {Math.max(totalPages, 1)}
                  </span>
                  <button
                    type="button"
                    disabled={page >= totalPages - 1}
                    onClick={() => setPage((p) => p + 1)}
                    className="rounded-md border border-border px-3 py-1 text-xs font-medium transition-colors hover:bg-muted disabled:opacity-40 disabled:cursor-not-allowed"
                  >
                    Next
                  </button>
                </div>
              </div>
            </div>
          </>
        ) : null}
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: Key Regulatory Findings                                   */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Key Regulatory Findings
        </h2>
        <p className="mb-4 text-sm leading-relaxed text-muted-foreground">
          The top disrupted transcription factors are master regulators of
          hepatocyte identity and metabolism, directly linking GWAS risk variants
          to disease-relevant regulatory circuits identified by SCENIC+.
        </p>
        <div className="grid gap-4 sm:grid-cols-3">
          {TOP_TFS.map((tf) => (
            <div
              key={tf.tf}
              className="rounded-lg border border-border bg-card px-4 py-4 shadow-sm"
            >
              <div className="flex items-center justify-between">
                <Link
                  href={`/gene?symbol=${encodeURIComponent(tf.tf)}`}
                  className="font-mono text-lg font-bold text-primary hover:underline"
                >
                  {tf.tf}
                </Link>
                <Badge variant="destructive" className="text-[10px]">
                  {tf.variants} variants
                </Badge>
              </div>
              <p className="mt-2 text-xs leading-relaxed text-muted-foreground">
                {tf.description}
              </p>
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/atlas"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          &larr; Atlas Construction
        </Link>
        <Link
          href="/explore"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Browse DEGs in Explorer &rarr;
        </Link>
        <Link
          href="/translation"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Therapeutic Translation &rarr;
        </Link>
      </div>
    </PageContainer>
  );
}
