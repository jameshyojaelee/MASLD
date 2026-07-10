"use client";

/**
 * Landing mini-viz: the Kleiner fibrosis-stage strip F0 -> F4.
 *
 * No fetch — a clean five-segment strip colored by the fibrosis-stage palette
 * (F0 = control gray "no fibrosis", F1-F4 warm progression). Segment labels are
 * drawn below the strip in theme-token text (never colored). An arrowhead marks
 * the direction of progression.
 */

import { FIBROSIS_STAGE_ORDER, fibrosisStageColor } from "@/lib/palette";
import { useMeasure } from "./use-measure";

const GAP = 3;

export function StageStrip() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const ready = width > 0 && height > 0;

  return (
    <div ref={ref} className="h-full w-full">
      {ready &&
        (() => {
          const n = FIBROSIS_STAGE_ORDER.length;
          const labelH = 14;
          const stripH = Math.max(6, Math.min(28, height - labelH - 4));
          const stripY = (height - labelH - stripH) / 2;
          const segW = (width - GAP * (n - 1)) / n;
          const fs = 9;
          return (
            <svg
              width={width}
              height={height}
              className="block"
              role="img"
              aria-label="Fibrosis stage progression F0 to F4"
            >
              {FIBROSIS_STAGE_ORDER.map((stage, i) => {
                const x = i * (segW + GAP);
                const cx = x + segW / 2;
                return (
                  <g key={stage}>
                    <rect
                      x={x}
                      y={stripY}
                      width={segW}
                      height={stripH}
                      rx={2}
                      fill={fibrosisStageColor(stage)}
                    />
                    <text
                      x={cx}
                      y={stripY + stripH + labelH - 3}
                      textAnchor="middle"
                      fontSize={fs}
                      fill="var(--color-muted-foreground)"
                      style={{ fontFamily: "var(--font-mono)" }}
                    >
                      {stage}
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
