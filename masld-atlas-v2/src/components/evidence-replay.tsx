"use client";

import { useState, useEffect, useMemo, useCallback } from "react";
import { EVIDENCE_SOURCES } from "@/lib/colors";
import type { EvidenceStrengths } from "@/lib/types";

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface ExpressionData {
  dream_logfc: number;
  dream_padj: number;
}

interface CausalData {
  coloc_pp4_max?: number;
  twas_pval?: number;
  twas_z?: number;
}

interface ProteomicsData {
  logfc?: number;
  padj?: number;
  tstat?: number;
  dataset?: string;
}

interface GeneProfileSnapshot {
  expression?: ExpressionData;
  causal?: CausalData;
  proteomics?: ProteomicsData;
}

interface EvidenceReplayProps {
  evidence: EvidenceStrengths;
  symbol: string;
  profile?: GeneProfileSnapshot;
  onClose: () => void;
}

// ---------------------------------------------------------------------------
// Narration text generator
// ---------------------------------------------------------------------------

function buildNarration(
  stepIndex: number,
  evidence: EvidenceStrengths,
  profile?: GeneProfileSnapshot
): { headline: string; detail: string } {
  const key = EVIDENCE_SOURCES[stepIndex]?.key;
  const val = key ? evidence[key] : 0;

  switch (stepIndex) {
    case 0: {
      const logfc = profile?.expression?.dream_logfc;
      const padj = profile?.expression?.dream_padj;
      if (val === 0) {
        return {
          headline: "Human Bulk RNA-seq",
          detail: "No significant differential expression in the Integrated mega-analysis.",
        };
      }
      const dir = logfc != null && logfc >= 0 ? "upregulated" : "downregulated";
      const lfcStr = logfc != null ? `logFC = ${logfc >= 0 ? "+" : ""}${logfc.toFixed(3)}` : "";
      const padjStr =
        padj != null
          ? padj < 1e-10
            ? "padj < 1e-10"
            : `padj = ${padj.toExponential(1)}`
          : "";
      return {
        headline: "Human Bulk RNA-seq",
        detail: `${dir.charAt(0).toUpperCase() + dir.slice(1)} across 10 cohorts (1,444 samples). ${lfcStr}${lfcStr && padjStr ? ", " : ""}${padjStr}`,
      };
    }
    case 1: {
      const pp4 = profile?.causal?.coloc_pp4_max;
      if (val === 0) {
        return {
          headline: "Genetic Causal",
          detail: "No colocalization detected across 24 GWAS studies.",
        };
      }
      const pp4Str = pp4 != null ? `Best COLOC PP.H4 = ${pp4.toFixed(3)}` : "";
      const strength =
        (pp4 ?? 0) >= 0.8 ? "strong" : (pp4 ?? 0) >= 0.5 ? "moderate" : "nominal";
      return {
        headline: "Genetic Causal",
        detail: `${strength.charAt(0).toUpperCase() + strength.slice(1)} genetic colocalization with MASLD GWAS. ${pp4Str}`,
      };
    }
    case 2:
      return val === 0
        ? { headline: "Essentiality", detail: "Not classified as essential in liver cell lines (DepMap CHRONOS)." }
        : {
            headline: "Essentiality",
            detail: `Liver-relevant essentiality signal detected (DepMap CHRONOS score: ${(val * 100).toFixed(0)}th percentile).`,
          };
    case 3:
      return val === 0
        ? { headline: "Epigenomic", detail: "No epigenomic regulation detected (SCENIC+ / scATAC-seq)." }
        : {
            headline: "Epigenomic",
            detail: `Epigenomic regulation confirmed. Active in SCENIC+ regulons and/or scATAC differential peaks. Strength: ${(val * 100).toFixed(0)}%.`,
          };
    case 4:
      return val === 0
        ? { headline: "Spatial", detail: "Not spatially variable in Visium liver sections." }
        : {
            headline: "Spatial",
            detail: `Spatially variable gene in Visium liver sections (Moran's I signal). Spatial score: ${(val * 100).toFixed(0)}%.`,
          };
    case 5:
      return val === 0
        ? { headline: "Single-Cell", detail: "Not differentially expressed in pseudobulk single-cell analysis." }
        : {
            headline: "Single-Cell",
            detail: `Significant in pseudobulk DE across cell types. Cell-type resolution score: ${(val * 100).toFixed(0)}%.`,
          };
    case 6:
      return val === 0
        ? { headline: "Mouse Concordance", detail: "No concordant signal across mouse diet models." }
        : {
            headline: "Mouse Concordance",
            detail: `Cross-species concordance confirmed across 5 mouse diet models. Concordance score: ${(val * 100).toFixed(0)}%.`,
          };
    case 7: {
      const lfc = profile?.proteomics?.logfc;
      const padj = profile?.proteomics?.padj;
      if (val === 0) {
        return {
          headline: "Proteomics",
          detail: "No proteomics signal in Olink plasma or DIA-MS liver datasets.",
        };
      }
      const dir = lfc != null && lfc >= 0 ? "Upregulated" : "Downregulated";
      const lfcStr = lfc != null ? `logFC = ${lfc >= 0 ? "+" : ""}${lfc.toFixed(2)}` : "";
      const padjStr =
        padj != null
          ? padj < 1e-10
            ? "padj < 1e-10"
            : `padj = ${padj.toExponential(1)}`
          : "";
      return {
        headline: "Proteomics",
        detail: `${dir} in proteomics (${lfcStr}${lfcStr && padjStr ? ", " : ""}${padjStr}).`,
      };
    }
    default:
      return { headline: "", detail: "" };
  }
}

// ---------------------------------------------------------------------------
// Radar SVG (reuses polygon math from evidence-fingerprint)
// ---------------------------------------------------------------------------

function ReplayRadar({
  evidence,
  step,
}: {
  evidence: EvidenceStrengths;
  step: number;
}) {
  const size = 300;
  const cx = size / 2;
  const cy = size / 2;
  const radius = size / 2 - 48; // extra padding for axis labels
  const n = EVIDENCE_SOURCES.length;

  const gridPath = useMemo(
    () =>
      Array.from({ length: n }, (_, i) => {
        const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
        return `${cx + radius * Math.cos(angle)},${cy + radius * Math.sin(angle)}`;
      }).join(" "),
    [cx, cy, radius, n]
  );

  // Points for the current step (cumulative fill)
  const cumulativePoints = useMemo(() => {
    return EVIDENCE_SOURCES.map((source, i) => {
      const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
      const val = i < step ? (evidence[source.key] ?? 0) : 0;
      const r = radius * val;
      return {
        x: cx + r * Math.cos(angle),
        y: cy + r * Math.sin(angle),
        val,
        source,
        angle,
        outerX: cx + radius * Math.cos(angle),
        outerY: cy + radius * Math.sin(angle),
        index: i,
      };
    });
  }, [evidence, step, cx, cy, radius, n]);

  const polygonPath = cumulativePoints.map((p) => `${p.x},${p.y}`).join(" ");

  // Color: use last-revealed source color, or primary if none
  const lastRevealedIndex = step - 1;
  const activeColor =
    lastRevealedIndex >= 0
      ? EVIDENCE_SOURCES[lastRevealedIndex].color
      : "var(--color-primary)";

  // Axis label positions (pushed outward beyond the outer grid)
  const labelRadius = radius + 30;

  return (
    <svg
      width={size}
      height={size}
      viewBox={`0 0 ${size} ${size}`}
      aria-label={`Radar chart showing ${step} of ${n} modalities`}
    >
      {/* Concentric grid rings */}
      {[0.25, 0.5, 0.75, 1].map((frac) => (
        <polygon
          key={frac}
          points={Array.from({ length: n }, (_, i) => {
            const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
            const r = radius * frac;
            return `${cx + r * Math.cos(angle)},${cy + r * Math.sin(angle)}`;
          }).join(" ")}
          fill="none"
          stroke="currentColor"
          strokeWidth={0.5}
          opacity={frac === 1 ? 0.2 : 0.08}
        />
      ))}

      {/* Spoke lines */}
      {cumulativePoints.map((p) => (
        <line
          key={p.source.key + "-spoke"}
          x1={cx}
          y1={cy}
          x2={p.outerX}
          y2={p.outerY}
          stroke="currentColor"
          strokeWidth={0.5}
          opacity={0.12}
        />
      ))}

      {/* Outer grid (full) */}
      <polygon
        points={gridPath}
        fill="none"
        stroke="currentColor"
        strokeWidth={1}
        opacity={0.2}
      />

      {/* Filled evidence polygon */}
      <polygon
        points={polygonPath}
        fill={activeColor}
        fillOpacity={0.2}
        stroke={activeColor}
        strokeWidth={2}
        strokeOpacity={0.8}
        style={{ transition: "all 0.6s ease" }}
      />

      {/* Vertex dots for revealed sources */}
      {cumulativePoints.map(
        (p) =>
          p.val > 0 && (
            <circle
              key={p.source.key + "-dot"}
              cx={p.x}
              cy={p.y}
              r={4}
              fill={p.source.color}
              stroke="var(--color-background)"
              strokeWidth={1.5}
              style={{ transition: "all 0.6s ease" }}
            />
          )
      )}

      {/* Axis labels */}
      {cumulativePoints.map((p) => {
        const lx = cx + labelRadius * Math.cos(p.angle);
        const ly = cy + labelRadius * Math.sin(p.angle);
        const isRevealed = p.index < step;
        return (
          <text
            key={p.source.key + "-label"}
            x={lx}
            y={ly}
            textAnchor="middle"
            dominantBaseline="middle"
            fontSize={9}
            fontWeight={isRevealed ? "600" : "400"}
            fill={isRevealed ? p.source.color : "currentColor"}
            opacity={isRevealed ? 1 : 0.3}
            style={{ transition: "all 0.6s ease" }}
          >
            {p.source.short}
          </text>
        );
      })}
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Progress dots
// ---------------------------------------------------------------------------

function ProgressDots({
  total,
  current,
  onClick,
}: {
  total: number;
  current: number;
  onClick: (i: number) => void;
}) {
  return (
    <div className="flex items-center gap-2">
      {Array.from({ length: total }, (_, i) => (
        <button
          key={i}
          onClick={() => onClick(i + 1)}
          aria-label={`Jump to step ${i + 1}: ${EVIDENCE_SOURCES[i].short}`}
          className="group relative"
        >
          <div
            className="h-2.5 w-2.5 rounded-full transition-all duration-300"
            style={{
              backgroundColor:
                i < current
                  ? EVIDENCE_SOURCES[i].color
                  : "currentColor",
              opacity: i < current ? 0.9 : 0.2,
              transform: i === current - 1 ? "scale(1.4)" : "scale(1)",
            }}
          />
        </button>
      ))}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main EvidenceReplay component
// ---------------------------------------------------------------------------

const TOTAL_STEPS = EVIDENCE_SOURCES.length;
const STEP_INTERVAL_MS = 1600;

export function EvidenceReplay({
  evidence,
  symbol,
  profile,
  onClose,
}: EvidenceReplayProps) {
  const [step, setStep] = useState(0);
  const [playing, setPlaying] = useState(false);

  // Count active sources
  const activeSources = useMemo(
    () => EVIDENCE_SOURCES.filter((s) => evidence[s.key] > 0).length,
    [evidence]
  );

  // Auto-advance when playing
  useEffect(() => {
    if (!playing) return;
    if (step >= TOTAL_STEPS) {
      setPlaying(false);
      return;
    }
    const timer = setTimeout(() => {
      setStep((s) => s + 1);
    }, STEP_INTERVAL_MS);
    return () => clearTimeout(timer);
  }, [playing, step]);

  const handlePlay = useCallback(() => {
    if (step >= TOTAL_STEPS) {
      // Restart
      setStep(0);
      setPlaying(true);
    } else {
      setPlaying((p) => !p);
    }
  }, [step]);

  const handleStepForward = useCallback(() => {
    setPlaying(false);
    setStep((s) => Math.min(s + 1, TOTAL_STEPS));
  }, []);

  const handleStepBack = useCallback(() => {
    setPlaying(false);
    setStep((s) => Math.max(s - 1, 0));
  }, []);

  const handleDotClick = useCallback((targetStep: number) => {
    setPlaying(false);
    setStep(targetStep);
  }, []);

  // Close on Escape
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, [onClose]);

  // Narration for the current step (0 = intro, 1-7 = sources)
  const narration = useMemo(() => {
    if (step === 0) {
      return {
        headline: symbol,
        detail: `Press play to watch evidence accumulate across ${TOTAL_STEPS} independent sources.`,
      };
    }
    if (step > TOTAL_STEPS) {
      return {
        headline: "Evidence Summary",
        detail: `${activeSources} of ${TOTAL_STEPS} modalities are active for ${symbol}.`,
      };
    }
    return buildNarration(step - 1, evidence, profile);
  }, [step, evidence, profile, symbol, activeSources]);

  const isIntro = step === 0;
  const isDone = step >= TOTAL_STEPS;

  return (
    /* Backdrop */
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/70 backdrop-blur-sm"
      onClick={onClose}
      role="dialog"
      aria-modal="true"
      aria-label={`Evidence replay for ${symbol}`}
    >
      {/* Panel — stop propagation so clicks inside don't close */}
      <div
        className="relative flex w-full max-w-lg flex-col items-center gap-6 rounded-2xl border border-border bg-background px-8 py-8 shadow-2xl"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Close button */}
        <button
          onClick={onClose}
          className="absolute right-4 top-4 flex h-8 w-8 items-center justify-center rounded-full text-muted-foreground transition-colors hover:bg-muted hover:text-foreground"
          aria-label="Close replay"
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="round"
          >
            <line x1="2" y1="2" x2="14" y2="14" />
            <line x1="14" y1="2" x2="2" y2="14" />
          </svg>
        </button>

        {/* Title */}
        <div className="text-center">
          <p className="text-xs font-semibold uppercase tracking-widest text-muted-foreground">
            Evidence Build
          </p>
          <h2 className="mt-0.5 font-mono text-2xl font-bold tracking-tight">
            {symbol}
          </h2>
        </div>

        {/* Radar */}
        <div className="flex items-center justify-center">
          <ReplayRadar evidence={evidence} step={step} />
        </div>

        {/* Narration */}
        <div
          className="min-h-[64px] w-full rounded-lg border border-border bg-muted/30 px-4 py-3 text-center"
          key={step} /* re-mount for fade-in via CSS animation */
          style={{ animation: "fadeInUp 0.35s ease" }}
        >
          {!isIntro && step <= TOTAL_STEPS && (
            <div
              className="mb-0.5 inline-block rounded px-2 py-0.5 text-xs font-semibold"
              style={{
                backgroundColor:
                  step > 0 && step <= TOTAL_STEPS
                    ? EVIDENCE_SOURCES[step - 1].color
                    : "transparent",
                color: "var(--color-background)",
                opacity: 0.9,
              }}
            >
              S{step}
            </div>
          )}
          <p className="text-sm font-semibold text-foreground">{narration.headline}</p>
          <p className="mt-0.5 text-xs leading-relaxed text-muted-foreground">
            {narration.detail}
          </p>
        </div>

        {/* Progress dots */}
        <ProgressDots
          total={TOTAL_STEPS}
          current={step}
          onClick={handleDotClick}
        />

        {/* Controls */}
        <div className="flex items-center gap-3">
          {/* Step back */}
          <button
            onClick={handleStepBack}
            disabled={step === 0}
            className="flex h-9 w-9 items-center justify-center rounded-full border border-border text-muted-foreground transition-colors hover:bg-muted hover:text-foreground disabled:cursor-not-allowed disabled:opacity-30"
            aria-label="Step back"
          >
            <svg width="14" height="14" viewBox="0 0 14 14" fill="currentColor">
              <polygon points="10,2 4,7 10,12" />
            </svg>
          </button>

          {/* Play / Pause */}
          <button
            onClick={handlePlay}
            className="flex h-12 w-12 items-center justify-center rounded-full bg-foreground text-background shadow-md transition-transform hover:scale-105 active:scale-95"
            aria-label={playing ? "Pause" : isDone ? "Replay" : "Play"}
          >
            {playing ? (
              /* Pause icon */
              <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor">
                <rect x="3" y="2" width="4" height="12" rx="1" />
                <rect x="9" y="2" width="4" height="12" rx="1" />
              </svg>
            ) : isDone ? (
              /* Replay icon */
              <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor">
                <path d="M8 3a5 5 0 1 0 4.546 2.914l1.54-.487A7 7 0 1 1 8 1V3z" />
                <polygon points="8,0 12,3 8,6" />
              </svg>
            ) : (
              /* Play icon */
              <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor">
                <polygon points="4,2 13,8 4,14" />
              </svg>
            )}
          </button>

          {/* Step forward */}
          <button
            onClick={handleStepForward}
            disabled={step >= TOTAL_STEPS}
            className="flex h-9 w-9 items-center justify-center rounded-full border border-border text-muted-foreground transition-colors hover:bg-muted hover:text-foreground disabled:cursor-not-allowed disabled:opacity-30"
            aria-label="Step forward"
          >
            <svg width="14" height="14" viewBox="0 0 14 14" fill="currentColor">
              <polygon points="4,2 10,7 4,12" />
            </svg>
          </button>
        </div>

        {/* Summary footer (always visible) */}
        <p className="text-center text-xs text-muted-foreground/60">
          {isDone
            ? `${activeSources}/${TOTAL_STEPS} sources active`
            : `Step ${step} of ${TOTAL_STEPS}`}
          {" · "}Press{" "}
          <kbd className="rounded border border-border px-1 py-0.5 font-mono text-[10px]">
            Esc
          </kbd>{" "}
          to close
        </p>
      </div>

      {/* Inline keyframe animation */}
      <style>{`
        @keyframes fadeInUp {
          from { opacity: 0; transform: translateY(6px); }
          to   { opacity: 1; transform: translateY(0); }
        }
      `}</style>
    </div>
  );
}
