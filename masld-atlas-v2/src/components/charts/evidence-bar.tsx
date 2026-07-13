"use client";

/**
 * Compact modality-evidence bars (S1–S8) for the gene card.
 *
 * One horizontal bar per evidence source, length = strength in [0, 1], colored
 * by the source's modality hue (MODALITY_HEX). Missing/zero modalities render as
 * an empty track so the eight-source grid is always legible at a glance.
 * Pure SVG-free markup, self-sizing, no axes — a small inline instrument.
 *
 * The bars reveal on mount (width 0 → strength, staggered per modality) via
 * framer's lightweight `m` — the previous CSS `transition-[width]` never fired
 * because the width was already at its final value on first paint. Reduced
 * motion renders the final widths instantly. When `onModalityClick` is provided,
 * each row becomes a button that jumps to that modality's detail section.
 */

import { useMemo } from "react";
import { m, useReducedMotion } from "framer-motion";
import { MODALITIES, type ModalityKey } from "@/lib/palette";

export interface EvidenceBarProps {
  /** Strength per modality, values in [0, 1]. Missing keys render empty. */
  strengths: Partial<Record<ModalityKey, number>>;
  /** Domain max for the bar length (default 1). */
  max?: number;
  /** Show the numeric value at the right of each bar. */
  showValues?: boolean;
  /** Use the short modality label (rendered to the left of each bar). */
  showLabels?: boolean;
  valueFormat?: (v: number) => string;
  className?: string;
  /** When set, each row becomes a button jumping to that modality's section. */
  onModalityClick?: (key: ModalityKey) => void;
}

export function EvidenceBar({
  strengths,
  max = 1,
  showValues = true,
  showLabels = true,
  valueFormat = (v) => v.toFixed(2),
  className,
  onModalityClick,
}: EvidenceBarProps) {
  const reduce = useReducedMotion();
  const rows = useMemo(
    () =>
      MODALITIES.map((m) => ({
        key: m.key,
        label: m.label,
        hex: m.hex,
        value: strengths[m.key] ?? 0,
      })),
    [strengths]
  );

  const labelW = showLabels ? 92 : 0;
  const valueW = showValues ? 40 : 0;

  return (
    <div className={className}>
      <div className="flex flex-col gap-1.5">
        {rows.map((r, i) => {
          const frac = Math.max(0, Math.min(1, max === 0 ? 0 : r.value / max));
          const clickable = !!onModalityClick;

          const inner = (
            <>
              {showLabels && (
                <span
                  className="shrink-0 text-right text-[11px] text-muted-foreground"
                  style={{ width: labelW }}
                >
                  {r.label}
                </span>
              )}
              <div className="relative h-2.5 flex-1 overflow-hidden rounded-full bg-muted">
                <m.div
                  className="absolute inset-y-0 left-0 rounded-full"
                  style={{ backgroundColor: r.hex, opacity: r.value > 0 ? 1 : 0 }}
                  initial={reduce ? false : { width: "0%" }}
                  animate={{ width: `${frac * 100}%` }}
                  transition={{
                    duration: 0.5,
                    delay: reduce ? 0 : i * 0.05,
                    ease: [0.16, 1, 0.3, 1],
                  }}
                />
              </div>
              {showValues && (
                <span
                  className="shrink-0 text-right text-[11px] tabular-nums text-muted-foreground"
                  style={{ width: valueW }}
                >
                  {r.value > 0 ? valueFormat(r.value) : "—"}
                </span>
              )}
            </>
          );

          return clickable ? (
            <button
              key={r.key}
              type="button"
              onClick={() => onModalityClick!(r.key)}
              className="flex w-full items-center gap-2 rounded-sm px-1 py-0.5 text-left transition-colors hover:bg-muted/60"
              title={`Jump to ${r.label} evidence`}
            >
              {inner}
            </button>
          ) : (
            <div key={r.key} className="flex items-center gap-2">
              {inner}
            </div>
          );
        })}
      </div>
    </div>
  );
}
