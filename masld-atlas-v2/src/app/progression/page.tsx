"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";

// ---------------------------------------------------------------------------
// Types matching progression_journey.json
// ---------------------------------------------------------------------------

interface GeneEntry {
  symbol: string;
  logfc: number;
  padj: number;
}

interface SwitchGene {
  symbol: string;
  logfc_early: number;
  logfc_late: number;
}

interface StageData {
  stage: string;
  label: string;
  n_samples: number;
  n_degs_up: number;
  n_degs_down: number;
  top_genes_up: GeneEntry[];
  top_genes_down: GeneEntry[];
  s1_fraction: number;
  s2_fraction: number;
  key_pathways: string[];
}

interface F2Switch {
  total_late_degs: number;
  early_degs: number;
  new_in_late: number;
  new_inflammatory_genes: number;
  key_switch_genes: SwitchGene[];
}

interface ClassifierMetrics {
  f_ge_3_auroc: number;
  nas_ge_5_auroc: number;
  fibrosis_qwk: number;
}

interface ProgressionData {
  stages: StageData[];
  f2_switch: F2Switch;
  classifier_metrics: ClassifierMetrics;
}

// ---------------------------------------------------------------------------
// Stage metadata
// ---------------------------------------------------------------------------

const STAGE_META: Record<string, { color: string; bgClass: string; borderClass: string; badgeVariant: "default" | "secondary" | "destructive" | "outline" }> = {
  F0: { color: "#10b981", bgClass: "bg-emerald-50 dark:bg-emerald-950/30", borderClass: "border-emerald-200 dark:border-emerald-800", badgeVariant: "secondary" },
  F1: { color: "#3b82f6", bgClass: "bg-blue-50 dark:bg-blue-950/30", borderClass: "border-blue-200 dark:border-blue-800", badgeVariant: "secondary" },
  F2: { color: "#f59e0b", bgClass: "bg-amber-50 dark:bg-amber-950/30", borderClass: "border-amber-200 dark:border-amber-800", badgeVariant: "destructive" },
  F3: { color: "#f97316", bgClass: "bg-orange-50 dark:bg-orange-950/30", borderClass: "border-orange-200 dark:border-orange-800", badgeVariant: "destructive" },
  F4: { color: "#ef4444", bgClass: "bg-red-50 dark:bg-red-950/30", borderClass: "border-red-200 dark:border-red-800", badgeVariant: "destructive" },
};

const STAGE_LABELS: Record<string, string> = {
  F0: "Healthy / No Fibrosis",
  F1: "Mild Fibrosis",
  F2: "Moderate Fibrosis \u2014 THE SWITCH",
  F3: "Severe Fibrosis",
  F4: "Cirrhosis",
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function formatLogFC(val: number): string {
  return val >= 0 ? `+${val.toFixed(2)}` : val.toFixed(2);
}

function logfcColor(val: number): string {
  if (val > 0) return "text-red-500 dark:text-red-400";
  if (val < 0) return "text-blue-500 dark:text-blue-400";
  return "text-muted-foreground";
}

function formatPadj(val: number): string {
  if (val === 0) return "<1e-300";
  if (val < 1e-300) return "<1e-300";
  return val.toExponential(1);
}

function isNamedGene(symbol: string): boolean {
  return !symbol.startsWith("ENSG0000");
}

// ---------------------------------------------------------------------------
// SVG: Timeline Bar
// ---------------------------------------------------------------------------

function TimelineBar({ stages }: { stages: StageData[] }) {
  const w = 600;
  const h = 80;
  const padX = 50;
  const spacing = (w - 2 * padX) / (stages.length - 1);
  const baseY = 40;

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="mx-auto w-full max-w-2xl"
      role="img"
      aria-label="Fibrosis stage timeline from F0 to F4"
    >
      {/* Connecting line */}
      <line
        x1={padX}
        y1={baseY}
        x2={w - padX}
        y2={baseY}
        stroke="currentColor"
        strokeOpacity={0.2}
        strokeWidth={2}
      />
      {stages.map((stage, i) => {
        const cx = padX + i * spacing;
        const meta = STAGE_META[stage.stage] ?? STAGE_META.F0;
        const isF2 = stage.stage === "F2";
        const r = isF2 ? 18 : 12;

        return (
          <g key={stage.stage}>
            {/* Highlight glow for F2 */}
            {isF2 && (
              <circle
                cx={cx}
                cy={baseY}
                r={24}
                fill={meta.color}
                opacity={0.12}
              />
            )}
            <circle
              cx={cx}
              cy={baseY}
              r={r}
              fill={meta.color}
              opacity={isF2 ? 1 : 0.8}
            />
            <text
              x={cx}
              y={baseY + 1}
              textAnchor="middle"
              dominantBaseline="central"
              fontSize={isF2 ? 11 : 9}
              fontWeight={isF2 ? 700 : 600}
              fill="white"
            >
              {stage.stage}
            </text>
            {/* Label below */}
            <text
              x={cx}
              y={baseY + (isF2 ? 34 : 26)}
              textAnchor="middle"
              fontSize={8}
              fill="currentColor"
              opacity={0.6}
            >
              {stage.label}
            </text>
            {/* SWITCH label for F2 */}
            {isF2 && (
              <text
                x={cx}
                y={baseY - 28}
                textAnchor="middle"
                fontSize={9}
                fontWeight={700}
                fill={meta.color}
              >
                THE SWITCH
              </text>
            )}
          </g>
        );
      })}
    </svg>
  );
}

// ---------------------------------------------------------------------------
// SVG: NMF Subtype Stacked Area Chart
// ---------------------------------------------------------------------------

function NMFAreaChart({ stages }: { stages: StageData[] }) {
  const margin = { left: 50, right: 20, top: 20, bottom: 36 };
  const w = 600;
  const h = 220;
  const plotW = w - margin.left - margin.right;
  const plotH = h - margin.top - margin.bottom;
  const n = stages.length;

  const scaleX = (i: number) => margin.left + (i / (n - 1)) * plotW;
  const scaleY = (frac: number) => margin.top + (1 - frac) * plotH;

  // S2 area: bottom 0 to s2_fraction
  const s2Points = stages.map((s, i) => ({
    x: scaleX(i),
    yTop: scaleY(s.s2_fraction),
    yBot: scaleY(0),
  }));

  // S1 area: from s2_fraction to 1.0
  const s1Points = stages.map((s, i) => ({
    x: scaleX(i),
    yTop: scaleY(1),
    yBot: scaleY(s.s2_fraction),
  }));

  function areaPath(points: { x: number; yTop: number; yBot: number }[]) {
    const top = points.map((p, i) => `${i === 0 ? "M" : "L"} ${p.x} ${p.yTop}`).join(" ");
    const bot = [...points].reverse().map((p, i) => `${i === 0 ? "L" : "L"} ${p.x} ${p.yBot}`).join(" ");
    return `${top} ${bot} Z`;
  }

  // Find approximate crossover (where S2 >= 0.5)
  const crossIdx = stages.findIndex((s) => s.s2_fraction >= 0.5);

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="mx-auto w-full max-w-2xl"
      role="img"
      aria-label="NMF subtype dynamics showing S1 and S2 fractions across fibrosis stages"
    >
      {/* Grid lines */}
      {[0, 0.25, 0.5, 0.75, 1].map((frac) => (
        <line
          key={frac}
          x1={margin.left}
          y1={scaleY(frac)}
          x2={w - margin.right}
          y2={scaleY(frac)}
          stroke="currentColor"
          strokeOpacity={0.08}
          strokeWidth={1}
        />
      ))}

      {/* S1 area (blue) */}
      <path d={areaPath(s1Points)} fill="#3b82f6" fillOpacity={0.25} stroke="#3b82f6" strokeWidth={1.5} strokeOpacity={0.6} />

      {/* S2 area (orange/red) */}
      <path d={areaPath(s2Points)} fill="#f97316" fillOpacity={0.25} stroke="#f97316" strokeWidth={1.5} strokeOpacity={0.6} />

      {/* Dividing line (S2 fraction) */}
      {stages.map((s, i) => {
        if (i === 0) return null;
        const prev = stages[i - 1];
        return (
          <line
            key={`div-${i}`}
            x1={scaleX(i - 1)}
            y1={scaleY(prev.s2_fraction)}
            x2={scaleX(i)}
            y2={scaleY(s.s2_fraction)}
            stroke="currentColor"
            strokeOpacity={0.3}
            strokeWidth={1}
            strokeDasharray="4,2"
          />
        );
      })}

      {/* Data points on the S2 line */}
      {stages.map((s, i) => (
        <circle
          key={s.stage}
          cx={scaleX(i)}
          cy={scaleY(s.s2_fraction)}
          r={3.5}
          fill="white"
          stroke="#f97316"
          strokeWidth={2}
        />
      ))}

      {/* Crossover marker */}
      {crossIdx >= 0 && (
        <g>
          <line
            x1={scaleX(crossIdx)}
            y1={margin.top}
            x2={scaleX(crossIdx)}
            y2={margin.top + plotH}
            stroke="#f97316"
            strokeOpacity={0.3}
            strokeWidth={1}
            strokeDasharray="3,3"
          />
        </g>
      )}

      {/* X-axis labels */}
      {stages.map((s, i) => (
        <text
          key={`x-${s.stage}`}
          x={scaleX(i)}
          y={h - 8}
          textAnchor="middle"
          fontSize={11}
          fontWeight={600}
          fill="currentColor"
          opacity={0.7}
        >
          {s.stage}
        </text>
      ))}

      {/* Y-axis labels */}
      {[0, 25, 50, 75, 100].map((pct) => (
        <text
          key={`y-${pct}`}
          x={margin.left - 8}
          y={scaleY(pct / 100) + 3}
          textAnchor="end"
          fontSize={9}
          fill="currentColor"
          opacity={0.5}
          fontFamily="monospace"
        >
          {pct}%
        </text>
      ))}

      {/* Legend */}
      <rect x={margin.left + 8} y={margin.top + 4} width={10} height={10} rx={2} fill="#3b82f6" fillOpacity={0.5} />
      <text x={margin.left + 22} y={margin.top + 13} fontSize={10} fill="currentColor" opacity={0.7}>
        S1 (Metabolic)
      </text>
      <rect x={margin.left + 120} y={margin.top + 4} width={10} height={10} rx={2} fill="#f97316" fillOpacity={0.5} />
      <text x={margin.left + 134} y={margin.top + 13} fontSize={10} fill="currentColor" opacity={0.7}>
        S2 (Inflammatory)
      </text>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// SVG: Subtype Bar (inline stacked bar for each stage card)
// ---------------------------------------------------------------------------

function SubtypeBar({ s1, s2 }: { s1: number; s2: number }) {
  const s1Pct = Math.round(s1 * 100);
  const s2Pct = Math.round(s2 * 100);

  return (
    <div className="flex items-center gap-2">
      <div className="flex h-4 flex-1 overflow-hidden rounded-full">
        <div
          className="bg-blue-400/60 transition-all"
          style={{ width: `${s1Pct}%` }}
          title={`S1: ${s1Pct}%`}
        />
        <div
          className="bg-orange-400/60 transition-all"
          style={{ width: `${s2Pct}%` }}
          title={`S2: ${s2Pct}%`}
        />
      </div>
      <span className="w-20 text-right font-mono text-xs text-muted-foreground">
        S2: {s2Pct}%
      </span>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Gene list (top 5 genes for a stage)
// ---------------------------------------------------------------------------

function GeneList({ genes, direction }: { genes: GeneEntry[]; direction: "up" | "down" }) {
  const display = genes.slice(0, 5);
  if (display.length === 0) return <span className="text-xs text-muted-foreground">None at this stage</span>;

  return (
    <div className="space-y-0.5">
      {display.map((g) => (
        <div key={g.symbol} className="flex items-center justify-between gap-2 text-xs">
          {isNamedGene(g.symbol) ? (
            <Link
              href={`/gene/${encodeURIComponent(g.symbol)}`}
              className="font-mono font-semibold text-primary hover:underline"
            >
              {g.symbol}
            </Link>
          ) : (
            <span className="font-mono text-muted-foreground">{g.symbol}</span>
          )}
          <span className={`font-mono ${logfcColor(g.logfc)}`}>
            {formatLogFC(g.logfc)}
          </span>
        </div>
      ))}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Stage Card
// ---------------------------------------------------------------------------

function StageCard({ stage }: { stage: StageData }) {
  const meta = STAGE_META[stage.stage] ?? STAGE_META.F0;
  const label = STAGE_LABELS[stage.stage] ?? stage.label;
  const isF2 = stage.stage === "F2";
  const totalDegs = stage.n_degs_up + stage.n_degs_down;

  return (
    <div
      className={`rounded-lg border p-5 ${meta.bgClass} ${meta.borderClass} ${isF2 ? "ring-2 ring-amber-400/50 dark:ring-amber-500/30" : ""}`}
    >
      <div className="mb-3 flex items-start justify-between gap-3">
        <div className="flex items-center gap-2.5">
          <Badge variant={meta.badgeVariant} className="text-xs font-bold">
            {stage.stage}
          </Badge>
          <div>
            <h3 className="text-sm font-semibold">{label}</h3>
            <p className="text-xs text-muted-foreground">
              n = {stage.n_samples.toLocaleString()} samples
            </p>
          </div>
        </div>
        <div className="text-right">
          <p className="text-lg font-bold tabular-nums">{totalDegs.toLocaleString()}</p>
          <p className="text-xs text-muted-foreground">stage-unique DEGs</p>
        </div>
      </div>

      {/* DEG counts */}
      {totalDegs > 0 && (
        <div className="mb-3 flex gap-4 text-xs">
          <span className="text-red-500 dark:text-red-400">
            {stage.n_degs_up.toLocaleString()} up
          </span>
          <span className="text-blue-500 dark:text-blue-400">
            {stage.n_degs_down.toLocaleString()} down
          </span>
        </div>
      )}

      {/* Subtype bar */}
      <div className="mb-3">
        <p className="mb-1 text-xs font-medium text-muted-foreground">NMF Subtype Composition</p>
        <SubtypeBar s1={stage.s1_fraction} s2={stage.s2_fraction} />
      </div>

      {/* Gene lists */}
      {totalDegs > 0 && (
        <>
          <Separator className="my-3" />
          <div className="grid gap-4 sm:grid-cols-2">
            <div>
              <p className="mb-1.5 text-xs font-medium text-red-600 dark:text-red-400">
                Top Upregulated
              </p>
              <GeneList genes={stage.top_genes_up} direction="up" />
            </div>
            <div>
              <p className="mb-1.5 text-xs font-medium text-blue-600 dark:text-blue-400">
                Top Downregulated
              </p>
              <GeneList genes={stage.top_genes_down} direction="down" />
            </div>
          </div>
        </>
      )}

      {/* F2 special: no DEGs message */}
      {isF2 && totalDegs === 0 && (
        <>
          <Separator className="my-3" />
          <p className="text-sm text-muted-foreground italic">
            F2 is the inflection point. No stage-unique DEGs detected here because this
            stage marks the transition boundary between early (F0-F1) and late (F3-F4)
            disease programs.
          </p>
        </>
      )}

      {/* Pathways */}
      {stage.key_pathways.length > 0 && (
        <>
          <Separator className="my-3" />
          <p className="mb-1.5 text-xs font-medium text-muted-foreground">Key Pathways</p>
          <div className="flex flex-wrap gap-1">
            {stage.key_pathways.map((pw) => (
              <Badge key={pw} variant="outline" className="text-[10px]">
                {pw}
              </Badge>
            ))}
          </div>
        </>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// F2 Switch Highlight Section
// ---------------------------------------------------------------------------

function F2SwitchSection({ data }: { data: F2Switch }) {
  return (
    <div className="rounded-lg border-2 border-amber-400 bg-amber-50/80 p-6 dark:border-amber-600 dark:bg-amber-950/30">
      <div className="mb-4 flex items-center gap-3">
        <div className="flex h-10 w-10 items-center justify-center rounded-full bg-amber-500 text-white">
          <svg width="20" height="20" viewBox="0 0 20 20" fill="none">
            <path d="M10 2L12.5 8H17.5L13.5 12L15 18L10 14L5 18L6.5 12L2.5 8H7.5L10 2Z" fill="currentColor" />
          </svg>
        </div>
        <div>
          <h2 className="text-lg font-bold">The Metabolic-to-Inflammatory Switch</h2>
          <p className="text-sm text-muted-foreground">
            A discrete transition at F2 that defines irreversibility
          </p>
        </div>
      </div>

      {/* Key stats */}
      <div className="mb-4 grid grid-cols-2 gap-3 sm:grid-cols-4">
        <div className="rounded-md bg-white/60 p-3 dark:bg-white/5">
          <p className="text-2xl font-bold tabular-nums">{data.early_degs.toLocaleString()}</p>
          <p className="text-xs text-muted-foreground">Early DEGs (F0-F1)</p>
        </div>
        <div className="rounded-md bg-white/60 p-3 dark:bg-white/5">
          <p className="text-2xl font-bold tabular-nums">{data.total_late_degs.toLocaleString()}</p>
          <p className="text-xs text-muted-foreground">Late DEGs (F3-F4)</p>
        </div>
        <div className="rounded-md bg-white/60 p-3 dark:bg-white/5">
          <p className="text-2xl font-bold tabular-nums text-amber-600 dark:text-amber-400">
            {data.new_in_late.toLocaleString()}
          </p>
          <p className="text-xs text-muted-foreground">New in late stages</p>
        </div>
        <div className="rounded-md bg-white/60 p-3 dark:bg-white/5">
          <p className="text-2xl font-bold tabular-nums text-red-600 dark:text-red-400">
            {data.new_inflammatory_genes.toLocaleString()}
          </p>
          <p className="text-xs text-muted-foreground">New inflammatory genes</p>
        </div>
      </div>

      {/* Visual: Early vs Late */}
      <div className="mb-4 flex items-center gap-4">
        <div className="flex-1 rounded-md bg-blue-100 p-3 text-center dark:bg-blue-900/30">
          <p className="text-sm font-semibold text-blue-700 dark:text-blue-300">Early (F0-F1)</p>
          <p className="text-xs text-blue-600/70 dark:text-blue-400/70">Metabolic programs</p>
          <p className="mt-1 text-xs text-muted-foreground">
            Cholesterol homeostasis, lipid metabolism
          </p>
        </div>
        <div className="flex flex-col items-center">
          <svg width="40" height="24" viewBox="0 0 40 24" className="text-amber-500">
            <path d="M4 12H36M36 12L28 6M36 12L28 18" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round" fill="none" />
          </svg>
          <span className="text-[10px] font-bold text-amber-600 dark:text-amber-400">F2</span>
        </div>
        <div className="flex-1 rounded-md bg-red-100 p-3 text-center dark:bg-red-900/30">
          <p className="text-sm font-semibold text-red-700 dark:text-red-300">Late (F3-F4)</p>
          <p className="text-xs text-red-600/70 dark:text-red-400/70">Inflammatory programs</p>
          <p className="mt-1 text-xs text-muted-foreground">
            NF-kB signaling, EMT, ECM remodeling
          </p>
        </div>
      </div>

      {/* Switch genes */}
      <Separator className="my-3" />
      <p className="mb-2 text-xs font-medium text-muted-foreground">
        Key Switch Genes (largest early-to-late change)
      </p>
      <div className="grid gap-1.5 sm:grid-cols-2 lg:grid-cols-3">
        {data.key_switch_genes.slice(0, 9).map((g) => (
          <div key={g.symbol} className="flex items-center justify-between rounded bg-white/50 px-2 py-1 text-xs dark:bg-white/5">
            {isNamedGene(g.symbol) ? (
              <Link
                href={`/gene/${encodeURIComponent(g.symbol)}`}
                className="font-mono font-semibold text-primary hover:underline"
              >
                {g.symbol}
              </Link>
            ) : (
              <span className="font-mono text-muted-foreground">{g.symbol}</span>
            )}
            <span className="flex gap-1.5 font-mono">
              <span className="text-muted-foreground">
                {g.logfc_early === 0 ? "0" : formatLogFC(g.logfc_early)}
              </span>
              <span className="text-muted-foreground/50">&rarr;</span>
              <span className={logfcColor(g.logfc_late)}>
                {formatLogFC(g.logfc_late)}
              </span>
            </span>
          </div>
        ))}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Classifier Metrics Section
// ---------------------------------------------------------------------------

function ClassifierSection({ metrics }: { metrics: ClassifierMetrics }) {
  const cards = [
    {
      label: "F\u22653 Detection",
      sublabel: "AUROC",
      value: metrics.f_ge_3_auroc.toFixed(3),
      description: "Advanced fibrosis detection",
    },
    {
      label: "NAS\u22655 Detection",
      sublabel: "AUROC",
      value: metrics.nas_ge_5_auroc.toFixed(3),
      description: "NASH activity detection",
    },
    {
      label: "Fibrosis Staging",
      sublabel: "Weighted Kappa",
      value: metrics.fibrosis_qwk.toFixed(3),
      description: "Ordinal classification",
    },
  ];

  return (
    <div className="grid gap-4 sm:grid-cols-3">
      {cards.map((c) => (
        <div
          key={c.label}
          className="rounded-lg border border-border bg-card p-4 text-center"
        >
          <p className="text-xs font-medium text-muted-foreground">{c.label}</p>
          <p className="mt-1 text-3xl font-bold tabular-nums tracking-tight">{c.value}</p>
          <p className="text-xs text-muted-foreground">{c.sublabel}</p>
          <p className="mt-2 text-xs text-muted-foreground/70">{c.description}</p>
        </div>
      ))}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main Page
// ---------------------------------------------------------------------------

export default function ProgressionPage() {
  const [data, setData] = useState<ProgressionData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetch("/data/progression_journey.json")
      .then((res) => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
      })
      .then((json: ProgressionData) => {
        setData(json);
        setLoading(false);
      })
      .catch((err) => {
        setError(err.message);
        setLoading(false);
      });
  }, []);

  if (loading) {
    return (
      <div className="w-full px-6 py-10">
        <p className="text-muted-foreground">Loading progression data...</p>
      </div>
    );
  }

  if (error || !data) {
    return (
      <div className="w-full px-6 py-10">
        <p className="text-destructive">Failed to load data: {error ?? "unknown error"}</p>
      </div>
    );
  }

  return (
    <div className="w-full px-6 py-8">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Disease Progression
        </h1>
        <p className="mt-2 max-w-2xl text-muted-foreground">
          From steatosis to cirrhosis: a five-stage journey through MASLD
          fibrosis. Each stage is characterized by distinct transcriptomic
          programs, with a discrete metabolic-to-inflammatory switch at F2.
        </p>
      </div>

      {/* Timeline Bar */}
      <section className="mb-10">
        <TimelineBar stages={data.stages} />
      </section>

      {/* Stage Cards */}
      <section className="mb-12 space-y-4">
        <h2 className="text-xl font-semibold tracking-tight">Stage Profiles</h2>
        <p className="text-sm text-muted-foreground">
          One-vs-rest differential expression at each fibrosis stage. Gene counts
          reflect stage-unique DEGs (padj &lt; 0.05, |logFC| &gt; 0.2).
        </p>
        <div className="grid gap-4 lg:grid-cols-2">
          {data.stages.map((stage) => (
            <StageCard key={stage.stage} stage={stage} />
          ))}
        </div>
      </section>

      {/* F2 Switch Highlight */}
      <section className="mb-12">
        <F2SwitchSection data={data.f2_switch} />
      </section>

      {/* NMF Subtype Dynamics */}
      <section className="mb-12">
        <h2 className="mb-2 text-xl font-semibold tracking-tight">
          NMF Subtype Dynamics
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Two molecular subtypes (NMF k=2, cophenetic = 0.989) shift across
          fibrosis stages. S1 (metabolic) dominates early disease; S2
          (inflammatory) rises through progression, reaching parity at F4.
        </p>
        <div className="rounded-lg border border-border bg-card p-4">
          <NMFAreaChart stages={data.stages} />
        </div>
      </section>

      {/* Classifier Metrics */}
      <section className="mb-12">
        <h2 className="mb-2 text-xl font-semibold tracking-tight">
          Classifier Performance
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Machine learning classifiers trained on stage-specific transcriptomic
          signatures. Leave-one-cohort-out cross-validation ensures
          generalizability across datasets.
        </p>
        <ClassifierSection metrics={data.classifier_metrics} />
      </section>

      {/* Bottom CTA */}
      <section className="rounded-lg border border-border bg-muted/30 px-6 py-5 text-center">
        <p className="mb-3 text-sm text-muted-foreground">
          Explore individual genes and their stage trajectories in the Gene Explorer.
        </p>
        <Link href="/explore">
          <Button variant="default">
            Explore genes in the Gene Explorer
          </Button>
        </Link>
      </section>
    </div>
  );
}
