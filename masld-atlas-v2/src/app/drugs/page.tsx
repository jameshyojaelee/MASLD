"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface FunnelStage {
  stage: string;
  count: number;
  description: string;
}

interface TopCompound {
  name: string;
  target: string | null;
  moa: string | null;
  composite_score: number;
  score_reversal: number;
  score_network: number;
  score_dgidb: number;
}

interface ValidatedDrug {
  drug: string;
  target: string;
  stage: string;
  moa: string;
  support: string;
  is_deg: boolean;
  dream_logfc: number | null;
}

interface SexStats {
  female_biased: number;
  male_biased: number;
  balanced: number;
  total_compounds: number;
}

interface DrugPipelineData {
  funnel: FunnelStage[];
  top_compounds: TopCompound[];
  validated_drugs: ValidatedDrug[];
  sex_stats: SexStats;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function stageBadgeClass(stage: string): string {
  if (stage.toLowerCase().includes("fda approved")) {
    return "bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300";
  }
  if (stage.toLowerCase().includes("phase 3")) {
    return "bg-blue-100 text-blue-800 dark:bg-blue-900/40 dark:text-blue-300";
  }
  if (stage.toLowerCase().includes("phase 2")) {
    return "bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300";
  }
  return "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400";
}

function supportBadgeClass(support: string): string {
  switch (support) {
    case "Strong":
      return "bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300";
    case "Moderate":
      return "bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300";
    default:
      return "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400";
  }
}

function formatLogFC(val: number | null | undefined): string {
  if (val == null) return "\u2014";
  return val >= 0 ? `+${val.toFixed(3)}` : val.toFixed(3);
}

function logfcColor(val: number | null | undefined): string {
  if (val == null) return "text-muted-foreground";
  if (val > 0) return "text-red-500 dark:text-red-400";
  if (val < 0) return "text-blue-500 dark:text-blue-400";
  return "text-muted-foreground";
}

// ---------------------------------------------------------------------------
// Funnel SVG
// ---------------------------------------------------------------------------

const FUNNEL_COLORS = [
  { fill: "#93c5fd", stroke: "#3b82f6" }, // blue-300 / blue-500
  { fill: "#6ee7b7", stroke: "#10b981" }, // emerald-300 / emerald-500
  { fill: "#fcd34d", stroke: "#f59e0b" }, // amber-300 / amber-500
  { fill: "#a78bfa", stroke: "#8b5cf6" }, // violet-300 / violet-500
];

function FunnelChart({ stages }: { stages: FunnelStage[] }) {
  if (stages.length === 0) return null;

  const svgW = 700;
  const svgH = 240;
  const padX = 20;
  const padY = 30;
  const barH = 48;
  const gapY = 12;

  const maxCount = Math.max(...stages.map((s) => s.count));

  // compute widths: proportional to count, but minimum 80px for readability
  const availW = svgW - padX * 2;
  const widths = stages.map((s) => Math.max(80, (s.count / maxCount) * availW));

  return (
    <svg
      viewBox={`0 0 ${svgW} ${svgH}`}
      className="w-full max-w-3xl"
      role="img"
      aria-label="Drug pipeline funnel showing 4 stages of filtering"
    >
      {stages.map((stage, i) => {
        const y = padY + i * (barH + gapY);
        const w = widths[i];
        const x = (svgW - w) / 2;
        const color = FUNNEL_COLORS[i % FUNNEL_COLORS.length];

        // connector trapezoid between stages
        const connector =
          i < stages.length - 1 ? (
            <polygon
              key={`conn-${i}`}
              points={`${x},${y + barH} ${x + w},${y + barH} ${(svgW - widths[i + 1]) / 2 + widths[i + 1]},${y + barH + gapY} ${(svgW - widths[i + 1]) / 2},${y + barH + gapY}`}
              fill={color.fill}
              opacity={0.25}
            />
          ) : null;

        return (
          <g key={i}>
            {connector}
            <rect
              x={x}
              y={y}
              width={w}
              height={barH}
              rx={6}
              fill={color.fill}
              stroke={color.stroke}
              strokeWidth={1.5}
            />
            {/* Count (large) */}
            <text
              x={svgW / 2}
              y={y + 20}
              textAnchor="middle"
              className="fill-foreground text-[15px] font-bold"
            >
              {stage.count.toLocaleString()}
            </text>
            {/* Stage label */}
            <text
              x={svgW / 2}
              y={y + 36}
              textAnchor="middle"
              className="fill-muted-foreground text-[11px]"
            >
              {stage.stage}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Composite score bar (inline)
// ---------------------------------------------------------------------------

function ScoreBar({ value, max = 1 }: { value: number; max?: number }) {
  const pct = Math.min(100, (value / max) * 100);
  return (
    <div className="flex items-center gap-2">
      <div className="h-2 w-20 rounded-full bg-muted">
        <div
          className="h-2 rounded-full bg-primary"
          style={{ width: `${pct}%` }}
        />
      </div>
      <span className="font-mono text-xs">{value.toFixed(3)}</span>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Sex stratification bar
// ---------------------------------------------------------------------------

function SexBar({ stats }: { stats: SexStats }) {
  const total = stats.total_compounds;
  const femalePct = (stats.female_biased / total) * 100;
  const malePct = (stats.male_biased / total) * 100;
  const balancedPct = (stats.balanced / total) * 100;

  return (
    <div className="space-y-3">
      {/* stacked bar */}
      <div className="flex h-6 w-full overflow-hidden rounded-full">
        {femalePct > 0 && (
          <div
            className="flex items-center justify-center bg-pink-400 text-[10px] font-semibold text-white dark:bg-pink-500"
            style={{ width: `${femalePct}%` }}
          >
            {stats.female_biased}
          </div>
        )}
        {balancedPct > 0 && (
          <div
            className="flex items-center justify-center bg-zinc-300 text-[10px] font-semibold text-zinc-700 dark:bg-zinc-600 dark:text-zinc-200"
            style={{ width: `${Math.max(balancedPct, 3)}%` }}
          >
            {stats.balanced}
          </div>
        )}
        {malePct > 0 && (
          <div
            className="flex items-center justify-center bg-blue-400 text-[10px] font-semibold text-white dark:bg-blue-500"
            style={{ width: `${malePct}%` }}
          >
            {stats.male_biased}
          </div>
        )}
      </div>
      {/* legend */}
      <div className="flex gap-4 text-xs text-muted-foreground">
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full bg-pink-400" />
          Female-biased ({stats.female_biased})
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full bg-zinc-300 dark:bg-zinc-600" />
          Balanced ({stats.balanced})
        </span>
        <span className="flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-2.5 rounded-full bg-blue-400" />
          Male-biased ({stats.male_biased})
        </span>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Page component
// ---------------------------------------------------------------------------

export default function DrugsPage() {
  const [data, setData] = useState<DrugPipelineData | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    fetch("/data/drug_pipeline.json")
      .then((r) => r.json())
      .then((d: DrugPipelineData) => {
        setData(d);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, []);

  if (loading) {
    return (
      <div className="w-full px-6 py-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Drug Pipeline Explorer
        </h1>
        <p className="mt-6 text-muted-foreground">Loading drug data...</p>
      </div>
    );
  }

  if (!data) {
    return (
      <div className="w-full px-6 py-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Drug Pipeline Explorer
        </h1>
        <p className="mt-6 text-muted-foreground">
          Failed to load drug pipeline data.
        </p>
      </div>
    );
  }

  // Deduplicate validated drugs by drug name for card display
  // (some drugs like Lanifibranor appear multiple times for different targets)
  const drugCardMap = new Map<string, ValidatedDrug[]>();
  for (const d of data.validated_drugs) {
    const existing = drugCardMap.get(d.drug) ?? [];
    existing.push(d);
    drugCardMap.set(d.drug, existing);
  }
  const drugCards = Array.from(drugCardMap.entries());

  return (
    <div className="w-full px-6 py-8">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Drug Pipeline Explorer
        </h1>
        <p className="mt-2 text-muted-foreground">
          Translating transcriptomic disease signatures into drug candidates
          through LINCS L1000 reversal, network proximity, and clinical
          validation across 133 CGP compounds.
        </p>
      </div>

      {/* Section 1: Funnel */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Discovery Funnel
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Compounds are filtered through successive evidence layers, from
          signature reversal to clinical validation.
        </p>
        <div className="rounded-lg border border-border bg-muted/30 p-6">
          <FunnelChart stages={data.funnel} />
          {/* descriptions below funnel */}
          <div className="mt-4 grid grid-cols-1 gap-2 sm:grid-cols-2 lg:grid-cols-4">
            {data.funnel.map((s, i) => (
              <div key={i} className="text-xs text-muted-foreground">
                <span
                  className="mr-1.5 inline-block h-2.5 w-2.5 rounded-sm"
                  style={{
                    backgroundColor:
                      FUNNEL_COLORS[i % FUNNEL_COLORS.length].fill,
                  }}
                />
                {s.description}
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Section 2: Top Reversal Compounds */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Top Reversal Compounds
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Top 20 compounds ranked by composite score (reversal 35% +
          significance 20% + MR 20% + DGIdb 15% + network 10%).
        </p>
        <div className="overflow-x-auto rounded-lg border border-border">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-border bg-muted/50">
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  #
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Compound
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Target
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  MoA
                </th>
                <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                  Composite
                </th>
                <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                  Reversal
                </th>
                <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                  Network
                </th>
              </tr>
            </thead>
            <tbody>
              {data.top_compounds.map((c, i) => (
                <tr
                  key={i}
                  className="border-b border-border/50 transition-colors hover:bg-muted/30"
                >
                  <td className="px-3 py-1.5 font-mono text-xs text-muted-foreground">
                    {i + 1}
                  </td>
                  <td className="px-3 py-1.5 font-medium">{c.name}</td>
                  <td className="px-3 py-1.5">
                    {c.target ? (
                      <Link
                        href={`/gene/${encodeURIComponent(c.target)}`}
                        className="font-mono text-xs font-semibold text-primary hover:underline"
                      >
                        {c.target}
                      </Link>
                    ) : (
                      <span className="text-xs text-muted-foreground">
                        {"\u2014"}
                      </span>
                    )}
                  </td>
                  <td className="px-3 py-1.5 text-xs text-muted-foreground">
                    {c.moa ?? "\u2014"}
                  </td>
                  <td className="px-3 py-1.5">
                    <ScoreBar value={c.composite_score} />
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {c.score_reversal.toFixed(3)}
                  </td>
                  <td className="px-3 py-1.5 text-right font-mono text-xs">
                    {c.score_network > 0
                      ? c.score_network.toFixed(3)
                      : "\u2014"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      {/* Section 3: Clinically Validated Drugs */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Clinically Validated Drugs
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          Drugs currently in clinical trials or approved for MASLD/NASH, with
          atlas transcriptomic support level.
        </p>
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {drugCards.map(([drugName, entries]) => {
            // Use first entry for shared fields
            const primary = entries[0];
            return (
              <div
                key={drugName}
                className="flex flex-col gap-2 rounded-lg border border-border p-4 transition-colors hover:bg-muted/30"
              >
                <div className="flex items-start justify-between gap-2">
                  <h3 className="text-sm font-semibold leading-tight">
                    {drugName}
                  </h3>
                  <span
                    className={`inline-flex shrink-0 items-center rounded-full px-2 py-0.5 text-[10px] font-medium ${stageBadgeClass(primary.stage)}`}
                  >
                    {primary.stage}
                  </span>
                </div>
                <p className="text-xs text-muted-foreground">{primary.moa}</p>
                {/* Targets */}
                <div className="flex flex-wrap items-center gap-1.5">
                  <span className="text-xs text-muted-foreground">
                    Targets:
                  </span>
                  {entries.map((e, j) => (
                    <Link
                      key={j}
                      href={`/gene/${encodeURIComponent(e.target)}`}
                      className="font-mono text-xs font-semibold text-primary hover:underline"
                    >
                      {e.target}
                    </Link>
                  ))}
                </div>
                {/* Support + DEG */}
                <div className="mt-auto flex items-center gap-2 pt-1">
                  <span
                    className={`inline-flex items-center rounded-full px-2 py-0.5 text-[10px] font-medium ${supportBadgeClass(primary.support)}`}
                  >
                    {primary.support}
                  </span>
                  {primary.is_deg && (
                    <Badge variant="default" className="text-[10px]">
                      DEG
                    </Badge>
                  )}
                  {primary.dream_logfc != null && (
                    <span
                      className={`font-mono text-[10px] ${logfcColor(primary.dream_logfc)}`}
                    >
                      logFC {formatLogFC(primary.dream_logfc)}
                    </span>
                  )}
                </div>
              </div>
            );
          })}
        </div>
      </section>

      {/* Section 4: Sex Stratification */}
      <section className="mb-10">
        <h2 className="mb-1 text-xl font-semibold tracking-tight">
          Sex-Stratified Reversal
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          CGP reversal compounds classified by sex-biased efficacy. Of{" "}
          {data.sex_stats.total_compounds} compounds,{" "}
          {data.sex_stats.female_biased} ({((data.sex_stats.female_biased / data.sex_stats.total_compounds) * 100).toFixed(0)}%) show
          female-biased reversal signatures, consistent with the predominantly
          female transcriptomic MASLD signature.
        </p>
        <div className="rounded-lg border border-border p-4">
          <SexBar stats={data.sex_stats} />
        </div>
      </section>
    </div>
  );
}
