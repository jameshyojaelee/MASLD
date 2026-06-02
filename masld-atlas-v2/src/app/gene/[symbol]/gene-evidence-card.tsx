"use client";

import { useEffect, useState, useCallback } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { EvidenceReplay } from "@/components/evidence-replay";
import { useAppStore } from "@/lib/store";
import { EVIDENCE_SOURCES } from "@/lib/colors";
import type { EvidenceStrengths } from "@/lib/types";

// ---------------------------------------------------------------------------
// Types for the per-gene JSON profile
// ---------------------------------------------------------------------------

interface CohortDE {
  dataset: string;
  logfc: number;
  pval: number;
}

interface StageEntry {
  stage: string;
  logfc: number;
  padj: number;
}

interface NasEntry {
  nas: string;
  logfc: number;
  padj: number;
}

interface PerCohortContrasts {
  disease_vs_control?: CohortDE[];
  MASH_vs_control?: CohortDE[];
  MASL_vs_control?: CohortDE[];
  MASH_vs_MASL?: CohortDE[];
}

type ContrastKey = keyof PerCohortContrasts;

const CONTRAST_LABELS: Record<ContrastKey, string> = {
  disease_vs_control: "Disease vs Healthy control",
  MASH_vs_control: "MASH vs Healthy control",
  MASL_vs_control: "MASL vs Healthy control",
  MASH_vs_MASL: "MASH vs MASL",
};

interface ExpressionData {
  dream_logfc: number;
  dream_padj: number;
  dream_tstat: number;
  /** Legacy: single-contrast per-cohort dream DE */
  per_cohort_dream?: CohortDE[];
  /** New: per-cohort DE per contrast family (Phase 7C) */
  per_cohort?: PerCohortContrasts;
  /** Primary: stage vs F0 control (requires Phase 7C R pipeline output). */
  stage_trajectory_vs_control?: StageEntry[];
  /** F_i vs Healthy controls (05d_fibrosis_vs_healthy_dream.R). */
  stage_trajectory_vs_healthy?: StageEntry[];
  /** Legacy one-vs-rest staging output — fallback if vs_control missing. */
  stage_trajectory_vs_rest?: StageEntry[];
  /** NAS 0-8 trajectory (one-vs-rest). */
  nas_trajectory?: NasEntry[];
  /** NAS_i vs NAS_0 baseline (05c_nas_vs_baseline_dream.R). */
  nas_trajectory_vs_nas0?: NasEntry[];
  /** NAS_i vs Healthy (diagnosis_harmonized == Control). */
  nas_trajectory_vs_healthy?: NasEntry[];
}

type NasControlVariant = "vs_nas0" | "vs_healthy";
type FibControlVariant = "vs_f0" | "vs_healthy";

interface CausalData {
  coloc_pp4_max?: number;
  coloc_best_gwas?: string;
  coloc_n_gwas_05?: number;
  coloc_n_gwas_08?: number;
  coloc_by_gwas?: Record<string, number>;
  twas_z?: number;
  twas_pval?: number;
}

interface PseudobulkEntry {
  cell_type: string;
  logfc: number;
  padj: number;
}

interface CelltypeData {
  attribution_class?: string;
  pseudobulk_de?: PseudobulkEntry[];
}

interface SexSubtypeData {
  sex_class?: string;
  logfc_female?: number;
  logfc_male?: number;
  padj_female?: number;
  padj_male?: number;
  nmf_subtype?: string;
  nmf_direction?: string;
  nmf_logfc?: number;
}

interface DrugEntry {
  drug: string;
  stage?: string;
  moa?: string;
  support?: string;
}

interface LincsEntry {
  name: string;
  score: number;
  moa?: string;
}

interface TherapeuticData {
  dgidb_druggable?: boolean;
  drugs?: DrugEntry[];
  lincs_compounds?: LincsEntry[];
  progression_class?: string;
  progression_drug?: string;
}

interface ExternalLinks {
  genecards?: string;
  opentargets?: string;
  gtex?: string;
  pubmed?: string;
}

interface ProteomicsData {
  logfc?: number;
  padj?: number;
  tstat?: number;
  dataset?: string;
}

interface GeneProfile {
  symbol: string;
  ensembl_id: string;
  biotype: string;
  expression?: ExpressionData;
  causal?: CausalData;
  celltype?: CelltypeData;
  sex_subtype?: SexSubtypeData;
  therapeutic?: TherapeuticData;
  external?: ExternalLinks;
  proteomics?: ProteomicsData;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function formatPadj(val: number | null | undefined): string {
  if (val == null) return "\u2014";
  if (val === 0) return "<1e-300";
  if (val < 1e-300) return "<1e-300";
  return val.toExponential(1);
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

function directionArrow(val: number | null | undefined): string {
  if (val == null) return "";
  return val > 0 ? "\u2191" : val < 0 ? "\u2193" : "\u2192";
}

function pp4Color(val: number): string {
  if (val >= 0.8) return "bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300";
  if (val >= 0.5) return "bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300";
  return "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400";
}

function supportColor(support: string | undefined): string {
  switch (support) {
    case "Strong":
      return "bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300";
    case "Moderate":
      return "bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300";
    default:
      return "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400";
  }
}

function biotypeLabel(biotype: string): string {
  return biotype.replace(/_/g, " ");
}

/** Build evidence strengths from the gene_index entry (fetched alongside profile) */
function buildEvidenceFromProfile(profile: GeneProfile): EvidenceStrengths {
  // We build approximate evidence from what we know
  const e: EvidenceStrengths = {
    s1_human: 0,
    s2_genetic: 0,
    s3_essential: 0,
    s4_epigenomic: 0,
    s5_spatial: 0,
    s6_singlecell: 0,
    s7_mouse: 0,
    s8_proteomics: 0,
  };

  // S1: Human expression
  if (profile.expression?.dream_padj != null) {
    const padj = profile.expression.dream_padj;
    const absLfc = Math.abs(profile.expression.dream_logfc ?? 0);
    if (padj <= 0.05 && absLfc >= 0.2) e.s1_human = Math.min(1, 0.5 + absLfc);
    else if (padj <= 0.05) e.s1_human = 0.3;
    else e.s1_human = Math.min(0.15, absLfc * 0.5);
  }

  // S2: Genetic
  if (profile.causal) {
    const pp4 = profile.causal.coloc_pp4_max ?? 0;
    e.s2_genetic = Math.min(1, pp4);
    if (profile.causal.twas_pval != null && profile.causal.twas_pval < 0.05) {
      e.s2_genetic = Math.min(1, e.s2_genetic + 0.2);
    }
  }

  // S6: Single-cell
  if (profile.celltype?.pseudobulk_de?.length) {
    const nSig = profile.celltype.pseudobulk_de.filter((d) => d.padj < 0.05).length;
    e.s6_singlecell = Math.min(1, nSig * 0.25);
  }

  // S8: Proteomics
  if (profile.proteomics?.logfc != null) {
    const lfc = profile.proteomics.logfc;
    const padj = profile.proteomics.padj ?? 1;
    if (padj < 0.05) {
      e.s8_proteomics = Math.min(1, Math.abs(lfc) / 2);
    } else {
      e.s8_proteomics = Math.min(0.15, Math.abs(lfc) * 0.3);
    }
  }

  return e;
}

// ---------------------------------------------------------------------------
// SVG chart components
// ---------------------------------------------------------------------------

function ForestPlot({ data }: { data: CohortDE[] }) {
  const margin = { left: 100, right: 50, top: 4, bottom: 4 };
  const barH = 18;
  const gap = 2;
  const h = data.length * (barH + gap) + margin.top + margin.bottom;
  const w = 500;
  const plotW = w - margin.left - margin.right;

  const maxAbs = Math.max(1, ...data.map((d) => Math.abs(d.logfc)));
  const scale = (v: number) => (v / maxAbs) * (plotW / 2);
  const centerX = margin.left + plotW / 2;

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="w-full max-w-lg"
      role="img"
      aria-label="Forest plot of per-cohort log fold changes"
    >
      {/* Center line */}
      <line
        x1={centerX}
        y1={margin.top}
        x2={centerX}
        y2={h - margin.bottom}
        stroke="currentColor"
        strokeOpacity={0.2}
        strokeWidth={1}
        strokeDasharray="2,2"
      />
      {data.map((d, i) => {
        const y = margin.top + i * (barH + gap);
        const barWidth = Math.abs(scale(d.logfc));
        const x = d.logfc >= 0 ? centerX : centerX - barWidth;
        const isSignificant = d.pval < 0.05;
        const fillColor = d.logfc >= 0
          ? isSignificant ? "var(--color-red-500, #ef4444)" : "var(--color-red-300, #fca5a5)"
          : isSignificant ? "var(--color-blue-500, #3b82f6)" : "var(--color-blue-300, #93c5fd)";

        return (
          <g key={d.dataset}>
            {/* Dataset label */}
            <text
              x={margin.left - 6}
              y={y + barH / 2 + 1}
              textAnchor="end"
              fontSize={10}
              fill="currentColor"
              opacity={0.7}
              dominantBaseline="middle"
            >
              {d.dataset}
            </text>
            {/* Bar */}
            <rect
              x={x}
              y={y + 2}
              width={Math.max(1, barWidth)}
              height={barH - 4}
              rx={2}
              fill={fillColor}
              opacity={isSignificant ? 0.85 : 0.45}
            />
            {/* Value label */}
            <text
              x={d.logfc >= 0 ? centerX + barWidth + 4 : centerX - barWidth - 4}
              y={y + barH / 2 + 1}
              textAnchor={d.logfc >= 0 ? "start" : "end"}
              fontSize={9}
              fill="currentColor"
              opacity={0.6}
              dominantBaseline="middle"
              fontFamily="monospace"
            >
              {d.logfc >= 0 ? "+" : ""}{d.logfc.toFixed(2)}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

function StageTrajectoryChart({
  data,
  variant = "vs_control",
}: {
  data: StageEntry[];
  variant?: "vs_control" | "vs_rest" | "vs_healthy";
}) {
  const stripSuffix = (label: string) => {
    if (variant === "vs_rest") return label.replace("_vs_rest", "");
    if (variant === "vs_healthy") return label.replace(/_vs_Healthy$/i, "");
    return label.replace("_vs_F0", "");
  };
  const margin = { left: 40, right: 20, top: 16, bottom: 28 };
  const w = 440;
  const h = 160;
  const plotW = w - margin.left - margin.right;
  const plotH = h - margin.top - margin.bottom;

  if (data.length === 0) return null;

  const maxAbs = Math.max(0.1, ...data.map((d) => Math.abs(d.logfc)));
  const scaleX = (i: number) => margin.left + (i / Math.max(1, data.length - 1)) * plotW;
  const scaleY = (v: number) => margin.top + plotH / 2 - (v / maxAbs) * (plotH / 2);

  const linePath = data
    .map((d, i) => `${i === 0 ? "M" : "L"} ${scaleX(i)} ${scaleY(d.logfc)}`)
    .join(" ");

  const zeroY = scaleY(0);

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="w-full max-w-md"
      role="img"
      aria-label="Stage trajectory chart"
    >
      {/* Zero line */}
      <line
        x1={margin.left}
        y1={zeroY}
        x2={w - margin.right}
        y2={zeroY}
        stroke="currentColor"
        strokeOpacity={0.15}
        strokeWidth={1}
      />
      {/* Line path */}
      <path
        d={linePath}
        fill="none"
        stroke="var(--color-primary)"
        strokeWidth={2}
        strokeLinejoin="round"
      />
      {/* Data points */}
      {data.map((d, i) => {
        const significant = d.padj < 0.05;
        return (
          <g key={d.stage}>
            <circle
              cx={scaleX(i)}
              cy={scaleY(d.logfc)}
              r={significant ? 5 : 3.5}
              fill={significant ? "var(--color-primary)" : "var(--color-background)"}
              stroke="var(--color-primary)"
              strokeWidth={significant ? 2 : 1.5}
            />
            {/* Stage label */}
            <text
              x={scaleX(i)}
              y={h - 6}
              textAnchor="middle"
              fontSize={9}
              fill="currentColor"
              opacity={0.6}
            >
              {stripSuffix(d.stage)}
            </text>
            {/* Value */}
            <text
              x={scaleX(i)}
              y={scaleY(d.logfc) - 10}
              textAnchor="middle"
              fontSize={8}
              fill="currentColor"
              opacity={significant ? 0.8 : 0.4}
              fontFamily="monospace"
            >
              {d.logfc >= 0 ? "+" : ""}{d.logfc.toFixed(2)}
            </text>
          </g>
        );
      })}
      {/* Y-axis label */}
      <text
        x={10}
        y={margin.top + plotH / 2}
        textAnchor="middle"
        fontSize={9}
        fill="currentColor"
        opacity={0.5}
        transform={`rotate(-90, 10, ${margin.top + plotH / 2})`}
      >
        logFC
      </text>
    </svg>
  );
}

function NasTrajectoryChart({
  data,
  variant = "vs_rest",
}: {
  data: NasEntry[];
  variant?: "vs_rest" | "vs_nas0" | "vs_healthy";
}) {
  const margin = { left: 40, right: 20, top: 16, bottom: 28 };
  const w = 440;
  const h = 160;
  const plotW = w - margin.left - margin.right;
  const plotH = h - margin.top - margin.bottom;

  if (data.length === 0) return null;

  const labelOf = (s: string) =>
    s
      .replace(/_vs_rest$/, "")
      .replace(/_vs_NAS0$/, "")
      .replace(/_vs_Healthy$/i, "")
      .replace(/^NAS/, "NAS ");
  void variant; // reserved for future per-variant styling
  const maxAbs = Math.max(0.1, ...data.map((d) => Math.abs(d.logfc)));
  const scaleX = (i: number) =>
    margin.left + (i / Math.max(1, data.length - 1)) * plotW;
  const scaleY = (v: number) =>
    margin.top + plotH / 2 - (v / maxAbs) * (plotH / 2);

  const linePath = data
    .map((d, i) => `${i === 0 ? "M" : "L"} ${scaleX(i)} ${scaleY(d.logfc)}`)
    .join(" ");

  const zeroY = scaleY(0);

  return (
    <svg
      viewBox={`0 0 ${w} ${h}`}
      className="w-full max-w-md"
      role="img"
      aria-label="NAS trajectory chart"
    >
      <line
        x1={margin.left}
        y1={zeroY}
        x2={w - margin.right}
        y2={zeroY}
        stroke="currentColor"
        strokeOpacity={0.15}
        strokeWidth={1}
      />
      <path
        d={linePath}
        fill="none"
        stroke="var(--color-primary)"
        strokeWidth={2}
        strokeLinejoin="round"
      />
      {data.map((d, i) => {
        const significant = d.padj < 0.05;
        return (
          <g key={d.nas}>
            <circle
              cx={scaleX(i)}
              cy={scaleY(d.logfc)}
              r={significant ? 5 : 3.5}
              fill={significant ? "var(--color-primary)" : "var(--color-background)"}
              stroke="var(--color-primary)"
              strokeWidth={significant ? 2 : 1.5}
            />
            <text
              x={scaleX(i)}
              y={h - 6}
              textAnchor="middle"
              fontSize={9}
              fill="currentColor"
              opacity={0.6}
            >
              {labelOf(d.nas)}
            </text>
            <text
              x={scaleX(i)}
              y={scaleY(d.logfc) - 10}
              textAnchor="middle"
              fontSize={8}
              fill="currentColor"
              opacity={significant ? 0.8 : 0.4}
              fontFamily="monospace"
            >
              {d.logfc >= 0 ? "+" : ""}
              {d.logfc.toFixed(2)}
            </text>
          </g>
        );
      })}
      <text
        x={10}
        y={margin.top + plotH / 2}
        textAnchor="middle"
        fontSize={9}
        fill="currentColor"
        opacity={0.5}
        transform={`rotate(-90, 10, ${margin.top + plotH / 2})`}
      >
        logFC
      </text>
    </svg>
  );
}

function ColocBar({ gwas, pp4 }: { gwas: string; pp4: number }) {
  const barWidth = Math.max(2, pp4 * 100);
  return (
    <div className="flex items-center gap-2 py-0.5">
      <span className="w-52 truncate text-xs text-muted-foreground" title={gwas}>
        {gwas}
      </span>
      <div className="relative h-3.5 flex-1 rounded-full bg-muted/50">
        <div
          className="absolute inset-y-0 left-0 rounded-full"
          style={{
            width: `${barWidth}%`,
            backgroundColor: pp4 >= 0.8 ? "var(--color-emerald-500, #10b981)" : pp4 >= 0.5 ? "var(--color-amber-500, #f59e0b)" : "var(--color-zinc-400, #a1a1aa)",
            opacity: 0.7,
          }}
        />
      </div>
      <span className="w-14 text-right font-mono text-xs text-muted-foreground">
        {pp4.toFixed(4)}
      </span>
    </div>
  );
}

function SexComparisonBars({
  logfcFemale,
  logfcMale,
  padjFemale,
  padjMale,
}: {
  logfcFemale?: number;
  logfcMale?: number;
  padjFemale?: number;
  padjMale?: number;
}) {
  if (logfcFemale == null && logfcMale == null) return null;
  const maxAbs = Math.max(0.1, Math.abs(logfcFemale ?? 0), Math.abs(logfcMale ?? 0));

  const renderBar = (label: string, val: number | undefined, padj: number | undefined) => {
    if (val == null) return null;
    const pct = (Math.abs(val) / maxAbs) * 50;
    const isSignificant = padj != null && padj < 0.05;
    const isPos = val >= 0;
    return (
      <div className="flex items-center gap-2">
        <span className="w-16 text-right text-xs text-muted-foreground">{label}</span>
        <div className="relative h-5 flex-1">
          {/* Center line */}
          <div className="absolute inset-y-0 left-1/2 w-px bg-border" />
          <div
            className="absolute top-0.5 h-4 rounded"
            style={{
              left: isPos ? "50%" : `${50 - pct}%`,
              width: `${pct}%`,
              backgroundColor: isPos
                ? isSignificant ? "var(--color-red-500, #ef4444)" : "var(--color-red-300, #fca5a5)"
                : isSignificant ? "var(--color-blue-500, #3b82f6)" : "var(--color-blue-300, #93c5fd)",
              opacity: isSignificant ? 0.8 : 0.4,
            }}
          />
        </div>
        <span className={`w-20 text-right font-mono text-xs ${logfcColor(val)}`}>
          {formatLogFC(val)}
          {padj != null && (
            <span className="ml-1 text-muted-foreground">
              {isSignificant ? "*" : ""}
            </span>
          )}
        </span>
      </div>
    );
  };

  return (
    <div className="space-y-1">
      {renderBar("Female", logfcFemale, padjFemale)}
      {renderBar("Male", logfcMale, padjMale)}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Section components
// ---------------------------------------------------------------------------

function SectionCard({
  title,
  children,
}: {
  title: string;
  children: React.ReactNode;
}) {
  return (
    <section className="rounded-lg border border-border bg-background p-5">
      <h2 className="mb-3 text-sm font-semibold uppercase tracking-wider text-muted-foreground">
        {title}
      </h2>
      {children}
    </section>
  );
}

function EmptyState({ label }: { label: string }) {
  return (
    <p className="text-sm text-muted-foreground/60">
      No {label} data available for this gene.
    </p>
  );
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

export function GeneEvidenceCard({ symbol }: { symbol: string }) {
  const [profile, setProfile] = useState<GeneProfile | null>(null);
  const [indexEvidence, setIndexEvidence] = useState<EvidenceStrengths | null>(null);
  const [loading, setLoading] = useState(true);
  const [notFound, setNotFound] = useState(false);
  const [copied, setCopied] = useState(false);
  const [showReplay, setShowReplay] = useState(false);
  const [contrastChoice, setContrastChoice] = useState<ContrastKey>("disease_vs_control");
  const [nasControlVariant, setNasControlVariant] =
    useState<NasControlVariant>("vs_nas0");
  const [fibControlVariant, setFibControlVariant] =
    useState<FibControlVariant>("vs_f0");

  const { compareGenes, addCompareGene, removeCompareGene } = useAppStore();
  const isCompared = compareGenes.includes(symbol);

  useEffect(() => {
    setLoading(true);
    setNotFound(false);

    // Fetch gene profile
    fetch(`/data/genes/${encodeURIComponent(symbol)}.json`)
      .then((res) => {
        if (!res.ok) throw new Error("not found");
        return res.json();
      })
      .then((data: GeneProfile) => {
        setProfile(data);
        setLoading(false);
      })
      .catch(() => {
        setNotFound(true);
        setLoading(false);
      });

    // Fetch evidence from index (for the fingerprint)
    fetch("/data/gene_index.json")
      .then((res) => res.json())
      .then((index: Array<{ symbol: string; evidence?: Record<string, number> }>) => {
        const entry = index.find((g) => g.symbol === symbol);
        if (entry?.evidence) {
          const raw = entry.evidence;
          // gene_index.json uses the M1-M7 aligned scheme
          // (s2_mouse, s3_genetic, s4_essential, s5_epigenomic, s6_spatial,
          // s7_singlecell). Map back to the legacy EvidenceStrengths keys.
          const full: EvidenceStrengths = {
            s1_human: raw.s1_human ?? 0,
            s2_genetic: raw.s3_genetic ?? raw.s2_genetic ?? 0,
            s3_essential: raw.s4_essential ?? raw.s3_essential ?? 0,
            s4_epigenomic: raw.s5_epigenomic ?? raw.s4_epigenomic ?? 0,
            s5_spatial: raw.s6_spatial ?? raw.s5_spatial ?? 0,
            s6_singlecell: raw.s7_singlecell ?? raw.s6_singlecell ?? 0,
            s7_mouse: raw.s2_mouse ?? raw.s7_mouse ?? 0,
            s8_proteomics: raw.s8_proteomics ?? 0,
          };
          setIndexEvidence(full);
        }
      })
      .catch(() => {
        // Silently fail -- we can compute from profile
      });
  }, [symbol]);

  const handleCopyLink = useCallback(() => {
    if (typeof window !== "undefined") {
      navigator.clipboard.writeText(window.location.href).then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      });
    }
  }, []);

  // Evidence: prefer index (precomputed), fall back to profile-derived
  const evidence = indexEvidence ?? (profile ? buildEvidenceFromProfile(profile) : null);

  // ----- Loading state -----
  if (loading) {
    return (
      <div className="mx-auto max-w-4xl px-6 py-10">
        <div className="animate-pulse space-y-4">
          <div className="h-8 w-40 rounded bg-muted" />
          <div className="h-4 w-64 rounded bg-muted" />
          <div className="h-40 rounded bg-muted" />
        </div>
      </div>
    );
  }

  // ----- 404 state -----
  if (notFound || !profile) {
    return (
      <div className="mx-auto max-w-4xl px-6 py-10">
        <h1 className="font-mono text-3xl font-bold tracking-tight">{symbol}</h1>
        <p className="mt-4 text-muted-foreground">
          No gene profile found for <span className="font-mono font-semibold">{symbol}</span>.
          This gene may not be present in the MASLD atlas.
        </p>
        <Link href="/explore" className="mt-4 inline-block">
          <Button variant="outline" size="sm">
            Back to Explorer
          </Button>
        </Link>
      </div>
    );
  }

  const expr = profile.expression;
  const causal = profile.causal;
  const celltype = profile.celltype;
  const sexSub = profile.sex_subtype;
  const therapeutic = profile.therapeutic;
  const external = profile.external;

  return (
    <div className="mx-auto max-w-4xl px-6 py-8">
      {/* ================================================================ */}
      {/* HEADER                                                           */}
      {/* ================================================================ */}
      <div className="sticky top-0 z-10 -mx-6 mb-6 border-b border-border bg-background/95 px-6 pb-4 pt-4 backdrop-blur-sm">
        <div className="flex items-start justify-between gap-4">
          <div className="flex items-center gap-4">
            {/* Evidence fingerprint */}
            {evidence && (
              <EvidenceFingerprint
                evidence={evidence}
                size={80}
                showTooltip={true}
              />
            )}
            <div>
              <div className="flex items-center gap-3">
                <h1 className="font-mono text-3xl font-bold tracking-tight">
                  {profile.symbol}
                </h1>
                <Badge variant="secondary" className="text-[10px]">
                  {biotypeLabel(profile.biotype)}
                </Badge>
              </div>
              <p className="mt-0.5 text-sm text-muted-foreground">
                {profile.ensembl_id}
              </p>
            </div>
          </div>

          {/* Action buttons */}
          <div className="flex shrink-0 flex-wrap gap-2">
            {evidence && (
              <Button
                variant="outline"
                size="sm"
                onClick={() => setShowReplay(true)}
              >
                Watch Evidence Build
              </Button>
            )}
            <Button
              variant={isCompared ? "default" : "outline"}
              size="sm"
              onClick={() =>
                isCompared ? removeCompareGene(symbol) : addCompareGene(symbol)
              }
            >
              {isCompared ? "Remove from Compare" : "Add to Compare"}
            </Button>
            <Button variant="outline" size="sm" onClick={handleCopyLink}>
              {copied ? "Copied!" : "Copy Link"}
            </Button>
          </div>
        </div>
      </div>

      {/* ================================================================ */}
      {/* SECTIONS                                                         */}
      {/* ================================================================ */}
      <div className="space-y-6">
        {/* -------------------------------------------------------------- */}
        {/* Section 1: Expression Evidence                                  */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="Expression Evidence">
          {expr ? (
            <div className="space-y-5">
              {/* Integrated DEG summary (dream mega-analysis, 10 cohorts) */}
              <div className="flex flex-wrap items-center gap-3">
                <span className="text-xs font-semibold uppercase tracking-wide text-muted-foreground">
                  Integrated DEG
                </span>
                <div
                  className={`inline-flex items-center gap-1.5 rounded-md px-3 py-1.5 font-mono text-sm font-semibold ${
                    expr.dream_padj != null && expr.dream_padj <= 0.05
                      ? expr.dream_logfc >= 0
                        ? "bg-red-100 text-red-800 dark:bg-red-900/30 dark:text-red-300"
                        : "bg-blue-100 text-blue-800 dark:bg-blue-900/30 dark:text-blue-300"
                      : "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400"
                  }`}
                  title="Integrated mega-analysis across 10 cohorts (1,444 samples)"
                >
                  <span className="text-base">{directionArrow(expr.dream_logfc)}</span>
                  <span>logFC {formatLogFC(expr.dream_logfc)}</span>
                </div>
                <span className="text-xs text-muted-foreground">
                  padj = {formatPadj(expr.dream_padj)}
                </span>
                {expr.dream_tstat != null && (
                  <span className="text-xs text-muted-foreground">
                    t = {expr.dream_tstat.toFixed(2)}
                  </span>
                )}
              </div>

              {/* 2x2 Trajectory grid: fibrosis (row 1) + NAS (row 2) */}
              <div className="grid grid-cols-1 gap-4 min-[800px]:grid-cols-2">
                {/* Row 1 col 1: Fibrosis stage vs control (toggle: vs F0 | vs Healthy) */}
                {(() => {
                  const f0Data = expr.stage_trajectory_vs_control;
                  const healthyData = expr.stage_trajectory_vs_healthy;
                  const activeData =
                    fibControlVariant === "vs_f0" ? f0Data : healthyData;
                  const hasAny =
                    (f0Data && f0Data.length > 0) ||
                    (healthyData && healthyData.length > 0);
                  return (
                    <div className="rounded-md border border-border/40 p-3">
                      <div className="mb-2 flex items-center justify-between gap-2">
                        <h3 className="text-xs font-medium text-muted-foreground">
                          Fibrosis stage vs baseline
                        </h3>
                        <div
                          role="radiogroup"
                          aria-label="Fibrosis baseline"
                          className="inline-flex overflow-hidden rounded-full border border-border text-[10px]"
                        >
                          <button
                            type="button"
                            role="radio"
                            aria-checked={fibControlVariant === "vs_f0"}
                            onClick={() => setFibControlVariant("vs_f0")}
                            className={`px-2 py-0.5 ${
                              fibControlVariant === "vs_f0"
                                ? "bg-foreground text-background"
                                : "text-muted-foreground hover:bg-muted/50"
                            }`}
                          >
                            vs F0
                          </button>
                          <button
                            type="button"
                            role="radio"
                            aria-checked={fibControlVariant === "vs_healthy"}
                            onClick={() => setFibControlVariant("vs_healthy")}
                            className={`px-2 py-0.5 border-l border-border ${
                              fibControlVariant === "vs_healthy"
                                ? "bg-foreground text-background"
                                : "text-muted-foreground hover:bg-muted/50"
                            }`}
                          >
                            vs Healthy
                          </button>
                        </div>
                      </div>
                      {activeData && activeData.length > 0 ? (
                        <>
                          <StageTrajectoryChart
                            data={activeData}
                            variant={
                              fibControlVariant === "vs_f0"
                                ? "vs_control"
                                : "vs_healthy"
                            }
                          />
                          <p className="mt-1 text-[10px] text-muted-foreground/60">
                            Filled circles indicate padj &lt; 0.05.
                          </p>
                        </>
                      ) : (
                        <p className="text-[11px] text-muted-foreground/70">
                          {fibControlVariant === "vs_f0"
                            ? "No stage vs F0 data for this gene."
                            : "No stage vs Healthy data for this gene."}
                          {!hasAny && " Awaiting pipeline output."}
                        </p>
                      )}
                    </div>
                  );
                })()}

                {/* Row 1 col 2: Fibrosis one-vs-rest */}
                <div className="rounded-md border border-border/40 p-3">
                  <h3 className="mb-2 text-xs font-medium text-muted-foreground">
                    Fibrosis stage (one-vs-rest)
                  </h3>
                  {expr.stage_trajectory_vs_rest && expr.stage_trajectory_vs_rest.length > 0 ? (
                    <>
                      <StageTrajectoryChart
                        data={expr.stage_trajectory_vs_rest}
                        variant="vs_rest"
                      />
                      <p className="mt-1 text-[10px] text-muted-foreground/60">
                        Filled circles indicate padj &lt; 0.05.
                      </p>
                    </>
                  ) : (
                    <p className="text-[11px] text-muted-foreground/70">
                      No one-vs-rest fibrosis data for this gene.
                    </p>
                  )}
                </div>

                {/* Row 2 col 1: NAS vs control (toggle between NAS_0 and Healthy) */}
                {(() => {
                  const nas0Data = expr.nas_trajectory_vs_nas0;
                  const healthyData = expr.nas_trajectory_vs_healthy;
                  const activeData =
                    nasControlVariant === "vs_nas0" ? nas0Data : healthyData;
                  const hasAny =
                    (nas0Data && nas0Data.length > 0) ||
                    (healthyData && healthyData.length > 0);
                  return (
                    <div className="rounded-md border border-border/40 p-3">
                      <div className="mb-2 flex items-center justify-between gap-2">
                        <h3 className="text-xs font-medium text-muted-foreground">
                          NAS vs baseline
                        </h3>
                        <div
                          role="radiogroup"
                          aria-label="NAS baseline"
                          className="inline-flex overflow-hidden rounded-full border border-border text-[10px]"
                        >
                          <button
                            type="button"
                            role="radio"
                            aria-checked={nasControlVariant === "vs_nas0"}
                            onClick={() => setNasControlVariant("vs_nas0")}
                            className={`px-2 py-0.5 ${
                              nasControlVariant === "vs_nas0"
                                ? "bg-foreground text-background"
                                : "text-muted-foreground hover:bg-muted/50"
                            }`}
                          >
                            vs NAS 0
                          </button>
                          <button
                            type="button"
                            role="radio"
                            aria-checked={nasControlVariant === "vs_healthy"}
                            onClick={() => setNasControlVariant("vs_healthy")}
                            className={`px-2 py-0.5 border-l border-border ${
                              nasControlVariant === "vs_healthy"
                                ? "bg-foreground text-background"
                                : "text-muted-foreground hover:bg-muted/50"
                            }`}
                          >
                            vs Healthy
                          </button>
                        </div>
                      </div>
                      {activeData && activeData.length > 0 ? (
                        <>
                          <NasTrajectoryChart
                            data={activeData}
                            variant={nasControlVariant}
                          />
                          <p className="mt-1 text-[10px] text-muted-foreground/60">
                            Filled circles indicate padj &lt; 0.05.
                          </p>
                        </>
                      ) : (
                        <p className="text-[11px] text-muted-foreground/70">
                          {nasControlVariant === "vs_nas0"
                            ? "No NAS vs NAS0 data for this gene."
                            : "No NAS vs Healthy data for this gene."}
                          {!hasAny && " Awaiting pipeline output."}
                        </p>
                      )}
                    </div>
                  );
                })()}

                {/* Row 2 col 2: NAS one-vs-rest */}
                <div className="rounded-md border border-border/40 p-3">
                  <h3 className="mb-2 text-xs font-medium text-muted-foreground">
                    NAS score (one-vs-rest)
                  </h3>
                  {expr.nas_trajectory && expr.nas_trajectory.length > 0 ? (
                    <>
                      <NasTrajectoryChart
                        data={expr.nas_trajectory}
                        variant="vs_rest"
                      />
                      <p className="mt-1 text-[10px] text-muted-foreground/60">
                        Filled circles indicate padj &lt; 0.05.
                      </p>
                    </>
                  ) : (
                    <p className="text-[11px] text-muted-foreground/70">
                      No one-vs-rest NAS data for this gene.
                    </p>
                  )}
                </div>
              </div>

              {/* Per-cohort logFC with contrast selector (below grid) */}
              {(() => {
                const contrastMap = expr.per_cohort;
                const legacyPerCohort = expr.per_cohort_dream;
                const contrastData = contrastMap?.[contrastChoice];
                const availableContrasts = (Object.keys(CONTRAST_LABELS) as ContrastKey[])
                  .filter((k) => (contrastMap?.[k]?.length ?? 0) > 0);
                const showSelector = availableContrasts.length > 0;
                if (!showSelector && !legacyPerCohort) return null;
                const effectiveData =
                  contrastData && contrastData.length > 0
                    ? contrastData
                    : !showSelector && legacyPerCohort
                      ? legacyPerCohort
                      : [];
                return (
                  <div>
                    <div className="mb-2 flex items-center gap-3">
                      <h3 className="text-xs font-medium text-muted-foreground">
                        Per-Cohort logFC
                      </h3>
                      {showSelector && (
                        <label className="flex items-center gap-1.5 text-xs text-muted-foreground">
                          <span>Contrast:</span>
                          <select
                            value={contrastChoice}
                            onChange={(e) => setContrastChoice(e.target.value as ContrastKey)}
                            className="rounded border border-border bg-background px-1.5 py-0.5 text-xs"
                          >
                            {availableContrasts.map((k) => (
                              <option key={k} value={k}>
                                {CONTRAST_LABELS[k]}
                              </option>
                            ))}
                          </select>
                        </label>
                      )}
                    </div>
                    {effectiveData.length > 0 ? (
                      <ForestPlot data={effectiveData} />
                    ) : (
                      <p className="text-[11px] text-muted-foreground/70">
                        No valid cohorts for this contrast.
                      </p>
                    )}
                    {!showSelector && legacyPerCohort && (
                      <p className="mt-1 text-[10px] text-muted-foreground/60">
                        Awaiting per-contrast pipeline output; showing primary Integrated contrast.
                      </p>
                    )}
                  </div>
                );
              })()}
            </div>
          ) : (
            <EmptyState label="expression" />
          )}
        </SectionCard>

        {/* -------------------------------------------------------------- */}
        {/* Section 2: Genetic Causal Evidence                             */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="Genetic Causal Evidence">
          {causal ? (
            <div className="space-y-4">
              {/* Best COLOC */}
              {causal.coloc_pp4_max != null && (
                <div className="flex flex-wrap items-center gap-3">
                  <div
                    className={`inline-flex items-center gap-1.5 rounded-md px-3 py-1.5 font-mono text-sm font-semibold ${pp4Color(causal.coloc_pp4_max)}`}
                  >
                    PP.H4 = {causal.coloc_pp4_max.toFixed(4)}
                  </div>
                  {causal.coloc_best_gwas && (
                    <span className="text-xs text-muted-foreground">
                      Best: {causal.coloc_best_gwas}
                    </span>
                  )}
                  {causal.coloc_n_gwas_08 != null && causal.coloc_n_gwas_08 > 0 && (
                    <Badge variant="secondary" className="text-[10px]">
                      {causal.coloc_n_gwas_08} GWAS PP4&gt;0.8
                    </Badge>
                  )}
                </div>
              )}

              {/* TWAS */}
              {causal.twas_z != null && (
                <div className="flex items-center gap-3">
                  <span className="text-xs font-medium text-muted-foreground">TWAS:</span>
                  <span className="font-mono text-xs">
                    z = {causal.twas_z.toFixed(2)}
                  </span>
                  <span className="text-xs text-muted-foreground">
                    p = {causal.twas_pval != null ? formatPadj(causal.twas_pval) : "\u2014"}
                  </span>
                  {causal.twas_pval != null && causal.twas_pval < 0.05 && (
                    <Badge variant="default" className="text-[10px]">
                      Significant
                    </Badge>
                  )}
                </div>
              )}

              {/* COLOC by GWAS */}
              {causal.coloc_by_gwas && Object.keys(causal.coloc_by_gwas).length > 0 && (
                <div>
                  <h3 className="mb-2 text-xs font-medium text-muted-foreground">
                    COLOC PP.H4 by GWAS ({Object.keys(causal.coloc_by_gwas).length} studies)
                  </h3>
                  <div className="space-y-0.5">
                    {Object.entries(causal.coloc_by_gwas)
                      .sort(([, a], [, b]) => b - a)
                      .map(([gwas, pp4]) => (
                        <ColocBar key={gwas} gwas={gwas} pp4={pp4} />
                      ))}
                  </div>
                </div>
              )}
            </div>
          ) : (
            <EmptyState label="genetic causal" />
          )}
        </SectionCard>

        {/* -------------------------------------------------------------- */}
        {/* Section 3: Cell-Type Resolution                                */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="Cell-Type Resolution">
          {celltype ? (
            <div className="space-y-3">
              {/* Attribution class */}
              {celltype.attribution_class && (
                <div>
                  <Badge
                    variant="secondary"
                    className="text-xs"
                  >
                    {celltype.attribution_class.replace(/_/g, " ")}
                  </Badge>
                </div>
              )}

              {/* Pseudobulk DE table */}
              {celltype.pseudobulk_de && celltype.pseudobulk_de.length > 0 ? (
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b border-border">
                        <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">
                          Cell Type
                        </th>
                        <th className="px-3 py-1.5 text-right text-xs font-medium text-muted-foreground">
                          logFC
                        </th>
                        <th className="px-3 py-1.5 text-right text-xs font-medium text-muted-foreground">
                          padj
                        </th>
                      </tr>
                    </thead>
                    <tbody>
                      {celltype.pseudobulk_de.map((entry) => {
                        const sig = entry.padj < 0.05;
                        return (
                          <tr
                            key={entry.cell_type}
                            className={`border-b border-border/30 ${sig ? "font-medium" : ""}`}
                          >
                            <td className="px-3 py-1 text-xs">{entry.cell_type}</td>
                            <td
                              className={`px-3 py-1 text-right font-mono text-xs ${logfcColor(entry.logfc)}`}
                            >
                              {formatLogFC(entry.logfc)}
                            </td>
                            <td className="px-3 py-1 text-right font-mono text-xs text-muted-foreground">
                              {formatPadj(entry.padj)}
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              ) : (
                !celltype.attribution_class && <EmptyState label="cell-type DE" />
              )}
            </div>
          ) : (
            <EmptyState label="cell-type" />
          )}
        </SectionCard>

        {/* -------------------------------------------------------------- */}
        {/* Section 4: Sex & Subtype                                       */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="Sex & Subtype">
          {sexSub ? (
            <div className="space-y-4">
              {/* Sex class */}
              {sexSub.sex_class && (
                <div>
                  <Badge
                    variant={
                      sexSub.sex_class.toLowerCase().includes("female")
                        ? "default"
                        : sexSub.sex_class.toLowerCase().includes("male")
                          ? "secondary"
                          : "outline"
                    }
                    className="text-xs"
                  >
                    {sexSub.sex_class.replace(/_/g, " ")}
                  </Badge>
                </div>
              )}

              {/* Sex logFC comparison */}
              <SexComparisonBars
                logfcFemale={sexSub.logfc_female}
                logfcMale={sexSub.logfc_male}
                padjFemale={sexSub.padj_female}
                padjMale={sexSub.padj_male}
              />

              {/* NMF subtype */}
              {sexSub.nmf_subtype && (
                <div className="flex items-center gap-3">
                  <Badge variant="outline" className="text-xs">
                    {sexSub.nmf_subtype} marker
                  </Badge>
                  {sexSub.nmf_direction && (
                    <span className="text-xs text-muted-foreground">
                      Direction: {sexSub.nmf_direction}
                    </span>
                  )}
                  {sexSub.nmf_logfc != null && (
                    <span className={`font-mono text-xs ${logfcColor(sexSub.nmf_logfc)}`}>
                      logFC {formatLogFC(sexSub.nmf_logfc)}
                    </span>
                  )}
                </div>
              )}

              {/* If only logfc values exist with no sex_class or nmf */}
              {!sexSub.sex_class && !sexSub.nmf_subtype && sexSub.logfc_female == null && sexSub.logfc_male == null && (
                <EmptyState label="sex/subtype" />
              )}
            </div>
          ) : (
            <EmptyState label="sex/subtype" />
          )}
        </SectionCard>

        {/* -------------------------------------------------------------- */}
        {/* Section 5: Therapeutic Actionability                            */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="Therapeutic Actionability">
          {therapeutic ? (
            <div className="space-y-4">
              {/* Druggability + progression class */}
              <div className="flex flex-wrap items-center gap-2">
                {therapeutic.dgidb_druggable != null && (
                  <Badge
                    variant={therapeutic.dgidb_druggable ? "default" : "secondary"}
                    className={`text-xs ${therapeutic.dgidb_druggable ? "bg-emerald-600 text-white dark:bg-emerald-700" : ""}`}
                  >
                    {therapeutic.dgidb_druggable ? "DGIdb Druggable" : "Not in DGIdb"}
                  </Badge>
                )}
                {therapeutic.progression_class && (
                  <Badge variant="outline" className="text-xs">
                    {therapeutic.progression_class.replace(/_/g, " ")}
                  </Badge>
                )}
              </div>

              {/* Drug table */}
              {therapeutic.drugs && therapeutic.drugs.length > 0 && (
                <div>
                  <h3 className="mb-2 text-xs font-medium text-muted-foreground">
                    Known Drugs ({therapeutic.drugs.length})
                  </h3>
                  <div className="overflow-x-auto rounded-md border border-border/50">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="border-b border-border bg-muted/30">
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">Drug</th>
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">Stage</th>
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">MoA</th>
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">Support</th>
                        </tr>
                      </thead>
                      <tbody>
                        {therapeutic.drugs.map((drug, i) => (
                          <tr key={i} className="border-b border-border/30">
                            <td className="px-3 py-1.5 text-xs font-medium">{drug.drug}</td>
                            <td className="px-3 py-1.5 text-xs text-muted-foreground">{drug.stage ?? "\u2014"}</td>
                            <td className="px-3 py-1.5 text-xs text-muted-foreground">{drug.moa ?? "\u2014"}</td>
                            <td className="px-3 py-1.5">
                              {drug.support && (
                                <span className={`inline-block rounded px-1.5 py-0.5 text-[10px] font-medium ${supportColor(drug.support)}`}>
                                  {drug.support}
                                </span>
                              )}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              {/* LINCS compounds */}
              {therapeutic.lincs_compounds && therapeutic.lincs_compounds.length > 0 && (
                <div>
                  <h3 className="mb-2 text-xs font-medium text-muted-foreground">
                    LINCS Reversal Compounds ({therapeutic.lincs_compounds.length})
                  </h3>
                  <div className="overflow-x-auto rounded-md border border-border/50">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="border-b border-border bg-muted/30">
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">Compound</th>
                          <th className="px-3 py-1.5 text-right text-xs font-medium text-muted-foreground">Score</th>
                          <th className="px-3 py-1.5 text-left text-xs font-medium text-muted-foreground">MoA</th>
                        </tr>
                      </thead>
                      <tbody>
                        {therapeutic.lincs_compounds.slice(0, 10).map((cmpd, i) => (
                          <tr key={i} className="border-b border-border/30">
                            <td className="px-3 py-1.5 text-xs font-medium">{cmpd.name}</td>
                            <td className="px-3 py-1.5 text-right font-mono text-xs text-muted-foreground">
                              {cmpd.score.toFixed(3)}
                            </td>
                            <td className="px-3 py-1.5 text-xs text-muted-foreground">{cmpd.moa ?? "\u2014"}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                    {therapeutic.lincs_compounds.length > 10 && (
                      <p className="px-3 py-1.5 text-[10px] text-muted-foreground">
                        Showing top 10 of {therapeutic.lincs_compounds.length} compounds.
                      </p>
                    )}
                  </div>
                </div>
              )}

              {/* Fallback: no drugs and no LINCS */}
              {(!therapeutic.drugs || therapeutic.drugs.length === 0) &&
                (!therapeutic.lincs_compounds || therapeutic.lincs_compounds.length === 0) &&
                !therapeutic.dgidb_druggable &&
                !therapeutic.progression_class && (
                <EmptyState label="therapeutic" />
              )}
            </div>
          ) : (
            <EmptyState label="therapeutic" />
          )}
        </SectionCard>

        {/* -------------------------------------------------------------- */}
        {/* Section 6: External Links                                      */}
        {/* -------------------------------------------------------------- */}
        <SectionCard title="External Resources">
          {external ? (
            <div className="flex flex-wrap gap-2">
              {external.genecards && (
                <a href={external.genecards} target="_blank" rel="noopener noreferrer">
                  <Button variant="outline" size="sm">
                    GeneCards
                  </Button>
                </a>
              )}
              {external.opentargets && (
                <a href={external.opentargets} target="_blank" rel="noopener noreferrer">
                  <Button variant="outline" size="sm">
                    Open Targets
                  </Button>
                </a>
              )}
              {external.gtex && (
                <a href={external.gtex} target="_blank" rel="noopener noreferrer">
                  <Button variant="outline" size="sm">
                    GTEx Portal
                  </Button>
                </a>
              )}
              {external.pubmed && (
                <a href={external.pubmed} target="_blank" rel="noopener noreferrer">
                  <Button variant="outline" size="sm">
                    PubMed
                  </Button>
                </a>
              )}
            </div>
          ) : (
            <EmptyState label="external link" />
          )}
        </SectionCard>
      </div>

      {/* Back link */}
      <div className="mt-8">
        <Link href="/explore">
          <Button variant="ghost" size="sm">
            &larr; Back to Explorer
          </Button>
        </Link>
      </div>

      {/* Evidence Replay modal */}
      {showReplay && evidence && (
        <EvidenceReplay
          evidence={evidence}
          symbol={profile.symbol}
          profile={{
            expression: expr
              ? { dream_logfc: expr.dream_logfc, dream_padj: expr.dream_padj }
              : undefined,
            causal: causal
              ? {
                  coloc_pp4_max: causal.coloc_pp4_max,
                  twas_pval: causal.twas_pval,
                  twas_z: causal.twas_z,
                }
              : undefined,
          }}
          onClose={() => setShowReplay(false)}
        />
      )}
    </div>
  );
}
