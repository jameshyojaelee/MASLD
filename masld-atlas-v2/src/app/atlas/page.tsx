"use client";

import { useEffect, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { getGeneIndex } from "@/lib/search-index";
import type { GeneIndexEntry } from "@/lib/types";

// ---------------------------------------------------------------------------
// Hardcoded study metadata
// ---------------------------------------------------------------------------

interface Cohort {
  dataset: string;
  samples: number;
  condition: string;
  platform: string;
  reference: string;
}

const COHORTS: Cohort[] = [
  {
    dataset: "GSE135251",
    samples: 216,
    condition: "MASLD staging (NAS)",
    platform: "RNA-seq",
    reference: "Govaere 2020",
  },
  {
    dataset: "GSE130970",
    samples: 78,
    condition: "NASH vs Healthy",
    platform: "RNA-seq",
    reference: "Hoang 2019",
  },
  {
    dataset: "GSE126848",
    samples: 112,
    condition: "NAFL/NASH/Fibrosis",
    platform: "RNA-seq",
    reference: "Suppli 2019",
  },
  {
    dataset: "GSE167523",
    samples: 118,
    condition: "NAFLD progression",
    platform: "RNA-seq",
    reference: "Pantano 2021",
  },
  {
    dataset: "GSE162694",
    samples: 76,
    condition: "Steatosis/NASH",
    platform: "RNA-seq",
    reference: "Fang 2021",
  },
  {
    dataset: "GSE174478",
    samples: 64,
    condition: "NAFL/NASH",
    platform: "RNA-seq",
    reference: "Vvedenskaya 2021",
  },
  {
    dataset: "GSE193066",
    samples: 216,
    condition: "MASLD staging",
    platform: "RNA-seq",
    reference: "Govaere 2022",
  },
  {
    dataset: "GSE213621",
    samples: 256,
    condition: "Fibrosis staging",
    platform: "RNA-seq",
    reference: "Vali 2023",
  },
  {
    dataset: "GSE240729",
    samples: 108,
    condition: "Fibrosis",
    platform: "RNA-seq",
    reference: "Liu 2024",
  },
  {
    dataset: "PRJNA512027",
    samples: 200,
    condition: "NASH/Healthy",
    platform: "RNA-seq",
    reference: "Gerhard 2018",
  },
];

const TOTAL_SAMPLES = 1444;

// ---------------------------------------------------------------------------
// Pipeline steps
// ---------------------------------------------------------------------------

interface PipelineStep {
  label: string;
  description: string;
}

const PIPELINE_STEPS: PipelineStep[] = [
  {
    label: "Per-study DE",
    description: "limma-voom per cohort with study-specific contrasts",
  },
  {
    label: "QC Filtering",
    description: "PCA + library size filters (pass_technical)",
  },
  {
    label: "Count Integration",
    description: "Harmonized count matrix across 10 cohorts",
  },
  {
    label: "Integrated Mega-Analysis",
    description: "Mixed-model mega-analysis with repeated measures",
  },
  {
    label: "DEG Thresholding",
    description: "Standard thresholds (padj < 0.05, |logFC| > 0.3)",
  },
  {
    label: "Consensus DEGs",
    description: "5,484 DEGs with LOO-CV validation (88.1% mean recovery)",
  },
];

// ---------------------------------------------------------------------------
// Key result cards
// ---------------------------------------------------------------------------

interface ResultCard {
  value: string;
  label: string;
  sublabel: string;
}

const RESULT_CARDS: ResultCard[] = [
  {
    value: "33,943",
    label: "Genes tested",
    sublabel: "GENCODE v49 / GRCh38",
  },
  {
    value: "5,484",
    label: "DEGs",
    sublabel: "padj < 0.05, |logFC| > 0.3",
  },
  {
    value: "88.1%",
    label: "LOO-CV recovery",
    sublabel: "rho = 0.959 across 10 folds",
  },
  {
    value: "7,548",
    label: "Robust genes",
    sublabel: "Present in 8/8 LOO folds",
  },
];

// ---------------------------------------------------------------------------
// Volcano SVG component
// ---------------------------------------------------------------------------

const VOLCANO_W = 700;
const VOLCANO_H = 480;
const MARGIN = { top: 20, right: 20, bottom: 50, left: 55 };
const INNER_W = VOLCANO_W - MARGIN.left - MARGIN.right;
const INNER_H = VOLCANO_H - MARGIN.top - MARGIN.bottom;

// Threshold lines
const LFC_THRESH = 0.3;
const PADJ_THRESH = 0.05;

interface VolcanoPoint {
  symbol: string;
  x: number; // logFC
  y: number; // -log10(padj)
  color: "up" | "down" | "ns";
}

function buildVolcanoPoints(genes: GeneIndexEntry[]): VolcanoPoint[] {
  const points: VolcanoPoint[] = [];
  for (const g of genes) {
    if (g.bulk_logfc == null || g.bulk_padj == null) continue;
    const lfc = g.bulk_logfc;
    const padj = g.bulk_padj;
    const neglog = Math.min(-Math.log10(Math.max(padj, 1e-300)), 300);
    let color: "up" | "down" | "ns" = "ns";
    if (g.is_deg && lfc > 0) color = "up";
    else if (g.is_deg && lfc < 0) color = "down";
    points.push({ symbol: g.symbol, x: lfc, y: neglog, color });
  }
  return points;
}

function scaleX(val: number, xMin: number, xMax: number): number {
  return ((val - xMin) / (xMax - xMin)) * INNER_W;
}

function scaleY(val: number, yMax: number): number {
  return INNER_H - (val / yMax) * INNER_H;
}

interface VolcanoPlotProps {
  genes: GeneIndexEntry[];
}

function VolcanoPlot({ genes }: VolcanoPlotProps) {
  const router = useRouter();
  const [tooltip, setTooltip] = useState<{
    symbol: string;
    lfc: number;
    padj: number;
    svgX: number;
    svgY: number;
  } | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  const points = buildVolcanoPoints(genes);

  if (points.length === 0) {
    return (
      <div className="flex h-[480px] items-center justify-center rounded-lg border border-border bg-muted/30 text-sm text-muted-foreground">
        Loading volcano plot...
      </div>
    );
  }

  const allX = points.map((p) => p.x);
  const allY = points.map((p) => p.y);
  const xMin = Math.min(...allX);
  const xMax = Math.max(...allX);
  const yMax = Math.max(...allY) || 10;

  // Threshold line positions
  const xThreshPos = scaleX(LFC_THRESH, xMin, xMax);
  const xThreshNeg = scaleX(-LFC_THRESH, xMin, xMax);
  const yThreshLine = scaleY(-Math.log10(PADJ_THRESH), yMax);

  // X axis tick values
  const xRange = xMax - xMin;
  const xTickStep = xRange > 4 ? 2 : 1;
  const xTicks: number[] = [];
  for (
    let t = Math.ceil(xMin / xTickStep) * xTickStep;
    t <= xMax;
    t += xTickStep
  ) {
    xTicks.push(parseFloat(t.toFixed(1)));
  }

  // Y axis tick values
  const yTickMax = Math.floor(yMax / 50) * 50;
  const yTicks: number[] = [0];
  for (let t = 50; t <= yTickMax; t += 50) yTicks.push(t);

  return (
    <div className="relative">
      <svg
        ref={svgRef}
        width={VOLCANO_W}
        height={VOLCANO_H}
        className="max-w-full overflow-visible"
        style={{ fontFamily: "inherit" }}
      >
        <g transform={`translate(${MARGIN.left},${MARGIN.top})`}>
          {/* Background */}
          <rect width={INNER_W} height={INNER_H} className="fill-muted/20" rx={4} />

          {/* Threshold lines */}
          <line
            x1={xThreshPos}
            x2={xThreshPos}
            y1={0}
            y2={INNER_H}
            stroke="hsl(var(--muted-foreground))"
            strokeWidth={0.8}
            strokeDasharray="4 3"
            opacity={0.5}
          />
          <line
            x1={xThreshNeg}
            x2={xThreshNeg}
            y1={0}
            y2={INNER_H}
            stroke="hsl(var(--muted-foreground))"
            strokeWidth={0.8}
            strokeDasharray="4 3"
            opacity={0.5}
          />
          <line
            x1={0}
            x2={INNER_W}
            y1={yThreshLine}
            y2={yThreshLine}
            stroke="hsl(var(--muted-foreground))"
            strokeWidth={0.8}
            strokeDasharray="4 3"
            opacity={0.5}
          />

          {/* Points — render non-significant first, then colored on top */}
          {points
            .filter((p) => p.color === "ns")
            .map((p, i) => {
              const cx = scaleX(p.x, xMin, xMax);
              const cy = scaleY(Math.min(p.y, yMax), yMax);
              return (
                <circle
                  key={`ns-${i}`}
                  cx={cx}
                  cy={cy}
                  r={1.5}
                  fill="hsl(var(--muted-foreground))"
                  opacity={0.45}
                />
              );
            })}
          {points
            .filter((p) => p.color !== "ns")
            .map((p, i) => {
              const cx = scaleX(p.x, xMin, xMax);
              const cy = scaleY(Math.min(p.y, yMax), yMax);
              const fill =
                p.color === "up"
                  ? "hsl(var(--destructive))"
                  : "hsl(210 100% 56%)";
              return (
                <circle
                  key={`sig-${i}`}
                  cx={cx}
                  cy={cy}
                  r={2.2}
                  fill={fill}
                  opacity={0.85}
                  style={{ cursor: "pointer" }}
                  onMouseEnter={() => {
                    const gene = genes.find((g) => g.symbol === p.symbol);
                    if (gene && gene.bulk_logfc != null && gene.bulk_padj != null) {
                      setTooltip({
                        symbol: p.symbol,
                        lfc: gene.bulk_logfc,
                        padj: gene.bulk_padj,
                        svgX: cx,
                        svgY: cy,
                      });
                    }
                  }}
                  onMouseLeave={() => setTooltip(null)}
                  onClick={() =>
                    router.push(`/gene/${encodeURIComponent(p.symbol)}`)
                  }
                />
              );
            })}

          {/* X axis */}
          <line x1={0} x2={INNER_W} y1={INNER_H} y2={INNER_H} stroke="hsl(var(--border))" strokeWidth={1} />
          {xTicks.map((t) => {
            const tx = scaleX(t, xMin, xMax);
            return (
              <g key={`xtick-${t}`}>
                <line x1={tx} x2={tx} y1={INNER_H} y2={INNER_H + 4} stroke="hsl(var(--border))" strokeWidth={1} />
                <text
                  x={tx}
                  y={INNER_H + 16}
                  textAnchor="middle"
                  fontSize={10}
                  fill="hsl(var(--muted-foreground))"
                >
                  {t.toFixed(1)}
                </text>
              </g>
            );
          })}
          <text
            x={INNER_W / 2}
            y={INNER_H + 38}
            textAnchor="middle"
            fontSize={11}
            fill="hsl(var(--muted-foreground))"
          >
            log\u2082FC (disease vs healthy)
          </text>

          {/* Y axis */}
          <line x1={0} x2={0} y1={0} y2={INNER_H} stroke="hsl(var(--border))" strokeWidth={1} />
          {yTicks.map((t) => {
            const ty = scaleY(t, yMax);
            return (
              <g key={`ytick-${t}`}>
                <line x1={-4} x2={0} y1={ty} y2={ty} stroke="hsl(var(--border))" strokeWidth={1} />
                <text
                  x={-8}
                  y={ty + 4}
                  textAnchor="end"
                  fontSize={10}
                  fill="hsl(var(--muted-foreground))"
                >
                  {t}
                </text>
              </g>
            );
          })}
          <text
            x={-INNER_H / 2}
            y={-40}
            textAnchor="middle"
            fontSize={11}
            fill="hsl(var(--muted-foreground))"
            transform="rotate(-90)"
          >
            -log\u2081\u2080(padj)
          </text>

          {/* Tooltip */}
          {tooltip && (() => {
            const tipW = 140;
            const tipH = 52;
            const tipX = tooltip.svgX + tipW > INNER_W ? tooltip.svgX - tipW - 6 : tooltip.svgX + 8;
            const tipY = tooltip.svgY - tipH < 0 ? tooltip.svgY + 6 : tooltip.svgY - tipH - 4;
            return (
              <g style={{ pointerEvents: "none" }}>
                <rect
                  x={tipX}
                  y={tipY}
                  width={tipW}
                  height={tipH}
                  rx={4}
                  fill="hsl(var(--popover))"
                  stroke="hsl(var(--border))"
                  strokeWidth={1}
                />
                <text x={tipX + 8} y={tipY + 17} fontSize={11} fontWeight={600} fill="hsl(var(--foreground))">
                  {tooltip.symbol}
                </text>
                <text x={tipX + 8} y={tipY + 31} fontSize={10} fill="hsl(var(--muted-foreground))">
                  logFC: {tooltip.lfc >= 0 ? "+" : ""}{tooltip.lfc.toFixed(3)}
                </text>
                <text x={tipX + 8} y={tipY + 44} fontSize={10} fill="hsl(var(--muted-foreground))">
                  padj: {tooltip.padj < 1e-300 ? "<1e-300" : tooltip.padj.toExponential(1)}
                </text>
              </g>
            );
          })()}
        </g>
      </svg>

      {/* Legend */}
      <div className="mt-2 flex flex-wrap gap-4 text-xs text-muted-foreground">
        <span className="flex items-center gap-1.5">
          <span
            className="inline-block h-2.5 w-2.5 rounded-full"
            style={{ background: "hsl(var(--destructive))", opacity: 0.75 }}
          />
          Up-regulated DEG
        </span>
        <span className="flex items-center gap-1.5">
          <span
            className="inline-block h-2.5 w-2.5 rounded-full"
            style={{ background: "hsl(210 100% 56%)", opacity: 0.75 }}
          />
          Down-regulated DEG
        </span>
        <span className="flex items-center gap-1.5">
          <span
            className="inline-block h-2.5 w-2.5 rounded-full"
            style={{ background: "hsl(var(--muted-foreground))", opacity: 0.3 }}
          />
          Non-significant
        </span>
        <span className="ml-auto italic">Click a colored dot to view gene</span>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function AtlasPage() {
  const [genes, setGenes] = useState<GeneIndexEntry[]>([]);
  const [volcanoLoading, setVolcanoLoading] = useState(true);

  useEffect(() => {
    getGeneIndex().then((data) => {
      setGenes(data);
      setVolcanoLoading(false);
    });
  }, []);

  return (
    <div className="mx-auto max-w-5xl px-6 py-10">
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">Atlas Construction</h1>
        <p className="mt-2 text-muted-foreground">
          10-cohort Integrated mega-analysis across 1,444 samples — 5,484 DEGs at
          padj&nbsp;&lt;&nbsp;0.05, |logFC|&nbsp;&gt;&nbsp;0.3, validated by
          leave-one-out cross-validation.
        </p>
      </div>

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Study Overview                                            */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Study Overview
        </h2>
        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                {["Dataset", "Samples", "Condition", "Platform", "Reference"].map(
                  (h) => (
                    <th
                      key={h}
                      className="px-4 py-2.5 text-left text-xs font-medium text-muted-foreground"
                    >
                      {h}
                    </th>
                  )
                )}
              </tr>
            </thead>
            <tbody>
              {COHORTS.map((c) => (
                <tr
                  key={c.dataset}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-4 py-2.5 font-mono text-xs font-semibold text-primary">
                    {c.dataset}
                  </td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs">
                    {c.samples}
                  </td>
                  <td className="px-4 py-2.5 text-xs">{c.condition}</td>
                  <td className="px-4 py-2.5">
                    <Badge variant="outline" className="text-[10px]">
                      {c.platform}
                    </Badge>
                  </td>
                  <td className="px-4 py-2.5 text-xs text-muted-foreground">
                    {c.reference}
                  </td>
                </tr>
              ))}
              {/* Total row */}
              <tr className="bg-muted/30">
                <td className="px-4 py-2.5 text-xs font-semibold">
                  Total (QC-passing)
                </td>
                <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold">
                  {TOTAL_SAMPLES.toLocaleString()}
                </td>
                <td colSpan={3} className="px-4 py-2.5 text-xs text-muted-foreground">
                  10 independent cohorts
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 2: Analysis Pipeline                                         */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Analysis Pipeline
        </h2>
        {/* Horizontal step flow */}
        <div className="flex flex-wrap items-start gap-0">
          {PIPELINE_STEPS.map((step, idx) => (
            <div key={step.label} className="flex items-start">
              {/* Step card */}
              <div className="flex w-[148px] flex-col gap-1.5 rounded-lg border border-border bg-card px-3 py-3 shadow-sm">
                <span className="text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
                  Step {idx + 1}
                </span>
                <span className="text-sm font-semibold leading-tight">
                  {step.label}
                </span>
                <span className="text-[11px] leading-tight text-muted-foreground">
                  {step.description}
                </span>
              </div>
              {/* Arrow connector */}
              {idx < PIPELINE_STEPS.length - 1 && (
                <div className="flex h-[72px] items-center px-1 text-muted-foreground">
                  <span className="text-lg">&#8594;</span>
                </div>
              )}
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: Key Results                                               */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Key Results
        </h2>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-4">
          {RESULT_CARDS.map((card) => (
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
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 4: Volcano Plot                                              */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <div className="mb-3 flex items-center justify-between">
          <h2 className="text-xl font-semibold tracking-tight">
            Volcano Plot
          </h2>
          {!volcanoLoading && (
            <span className="text-xs text-muted-foreground">
              {genes.filter((g) => g.bulk_logfc != null && g.bulk_padj != null).length.toLocaleString()} genes plotted
            </span>
          )}
        </div>
        <div className="rounded-lg border border-border bg-card p-4 shadow-sm">
          {volcanoLoading ? (
            <div className="flex h-[480px] items-center justify-center text-sm text-muted-foreground">
              Loading gene data&hellip;
            </div>
          ) : (
            <VolcanoPlot genes={genes} />
          )}
        </div>
        <p className="mt-2 text-xs text-muted-foreground">
          Dashed lines mark |logFC|&nbsp;=&nbsp;0.3 and padj&nbsp;=&nbsp;0.05
          thresholds. Colored dots are significant DEGs. Click any colored dot
          to view the gene profile. Points at y&nbsp;=&nbsp;300 represent
          padj&nbsp;&lt;&nbsp;1e-300 (clamped for display).
        </p>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/explore"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Browse DEGs in Explorer &rarr;
        </Link>
        <Link
          href="/causal"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Causal Architecture &rarr;
        </Link>
      </div>
    </div>
  );
}
