"use client";

/**
 * Landing mini-viz: the drug-development funnel (preclinical -> clinical ->
 * approved).
 *
 * No fetch — counts come from `atlas-constants`. Tier widths scale by log of
 * the count so the 2 approved targets stay visible next to 1,504 preclinical
 * ones. A single measured SVG draws three stacked trapezoids (approved darkest)
 * that fill the tile height, each with a compact inline label — count over tier
 * name — sitting beside its band. Stage colors come from the dev-stage palette
 * (preclinical = control gray); label text stays theme-token neutral (never
 * data-colored).
 */

import {
  DRUGS_APPROVED,
  DRUGS_CLINICAL,
  DRUGS_PRECLINICAL,
  fmt,
} from "@/lib/atlas-constants";
import { devStageColor } from "@/lib/palette";
import { useMeasure } from "./use-measure";

interface Tier {
  name: string;
  value: number;
  color: string;
}

const TIERS: Tier[] = [
  { name: "Preclinical", value: DRUGS_PRECLINICAL, color: devStageColor("Preclinical") },
  { name: "Clinical", value: DRUGS_CLINICAL, color: devStageColor("Clinical") },
  { name: "Approved", value: DRUGS_APPROVED, color: devStageColor("Approved") },
];

/** Minimum width fraction so the narrowest tier stays visibly clickable. */
const MIN_FRAC = 0.16;
/** Share of the width the funnel graphic claims (the rest holds the labels). */
const FUNNEL_FRAC = 0.55;

export function DrugFunnel() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const ready = width > 0 && height > 0;

  const maxLog = Math.max(...TIERS.map((t) => Math.log10(t.value + 1)));
  const frac = TIERS.map((t) =>
    Math.max(MIN_FRAC, Math.log10(t.value + 1) / maxLog)
  );

  return (
    <div ref={ref} className="h-full w-full">
      {ready &&
        (() => {
          const funnelW = width * FUNNEL_FRAC;
          const cx = funnelW / 2;
          const track = funnelW - 4;
          const bandH = height / TIERS.length;
          // Boundary widths: top of each band, plus a flat bottom.
          const w = frac.map((f) => f * track);
          const level = [...w, w[w.length - 1]];
          const labelX = funnelW + 8;
          const countFs = Math.min(12, Math.max(9, bandH * 0.34));
          const nameFs = Math.min(9, Math.max(7, bandH * 0.24));
          return (
            <svg
              width={width}
              height={height}
              className="block"
              role="img"
              aria-label="Drug development funnel: preclinical, clinical, approved"
            >
              {TIERS.map((t, k) => {
                const yTop = k * bandH + 1;
                const yBot = (k + 1) * bandH - 1;
                const yMid = (yTop + yBot) / 2;
                const tw = level[k];
                const bw = level[k + 1];
                const pts = [
                  [cx - tw / 2, yTop],
                  [cx + tw / 2, yTop],
                  [cx + bw / 2, yBot],
                  [cx - bw / 2, yBot],
                ]
                  .map((p) => p.join(","))
                  .join(" ");
                return (
                  <g key={t.name}>
                    <polygon points={pts} fill={t.color} />
                    <text
                      x={labelX}
                      y={yMid - nameFs * 0.55}
                      dominantBaseline="central"
                      fontSize={countFs}
                      fill="var(--color-foreground)"
                      style={{ fontFamily: "var(--font-mono)" }}
                    >
                      {fmt(t.value)}
                    </text>
                    <text
                      x={labelX}
                      y={yMid + countFs * 0.55}
                      dominantBaseline="central"
                      fontSize={nameFs}
                      fill="var(--color-muted-foreground)"
                    >
                      {t.name}
                    </text>
                  </g>
                );
              })}
            </svg>
          );
        })()}
    </div>
  );
}
