"use client";

/**
 * Landing mini-viz: the drug-development funnel (preclinical -> clinical ->
 * approved).
 *
 * No fetch — counts come from `atlas-constants`. Tier widths scale by log of
 * the count so the 2 approved targets stay visible next to 1,504 preclinical
 * ones. An SVG funnel (three stacked trapezoids, approved darkest) sits beside
 * a small theme-token legend; stage colors come from the dev-stage palette
 * (preclinical = control gray).
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

const MIN_FRAC = 0.16;

export function DrugFunnel() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const ready = width > 0 && height > 0;

  const maxLog = Math.max(...TIERS.map((t) => Math.log10(t.value + 1)));
  const frac = TIERS.map((t) =>
    Math.max(MIN_FRAC, Math.log10(t.value + 1) / maxLog)
  );

  return (
    <div className="flex h-full w-full items-center gap-3">
      <div ref={ref} className="h-full flex-[3]">
        {ready && (
          <svg
            width={width}
            height={height}
            className="block"
            role="img"
            aria-label="Drug development funnel"
          >
            {(() => {
              const cx = width / 2;
              const track = width - 2;
              const bandH = height / TIERS.length;
              // Boundary widths: top of each band, plus a flat bottom.
              const w = frac.map((f) => f * track);
              const level = [...w, w[w.length - 1]];
              return TIERS.map((t, k) => {
                const yTop = k * bandH + 1;
                const yBot = (k + 1) * bandH - 1;
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
                return <polygon key={t.name} points={pts} fill={t.color} />;
              });
            })()}
          </svg>
        )}
      </div>
      <ul className="flex flex-[2] flex-col justify-center gap-1.5">
        {TIERS.map((t) => (
          <li key={t.name} className="flex items-center gap-1.5 text-[10px] leading-tight">
            <span
              aria-hidden
              className="size-2 shrink-0 rounded-[2px]"
              style={{ backgroundColor: t.color }}
            />
            <span className="text-muted-foreground">{t.name}</span>
            <span className="ml-auto font-numeric text-foreground">{fmt(t.value)}</span>
          </li>
        ))}
      </ul>
    </div>
  );
}
