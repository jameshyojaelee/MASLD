"use client";

/**
 * Landing mini-viz: a tiny volcano scatter of the pooled bulk DEG contrast.
 *
 * Self-contained and fail-soft — fetches its own `volcano_preview.json`
 * (`[{symbol,x,y,sig}]`; x = logFC on a symmetric domain, y = -log10 padj) and
 * renders a skeleton while loading, nothing on failure. Marks follow the
 * project color contract: significant points are diverging red (up) / blue
 * (down) at full opacity; non-significant points are control gray at reduced
 * opacity. Text/gridlines use theme tokens (never a data hue).
 */

import { useEffect, useState } from "react";
import { dataUrl } from "@/lib/data-base";
import {
  DIVERGING,
  CONTROL,
  SIG_OPACITY,
  NONSIG_OPACITY,
} from "@/lib/palette";
import { useMeasure } from "./use-measure";

interface VolcanoPoint {
  symbol: string;
  x: number;
  y: number;
  sig: boolean;
}

const PAD = 6;

export function MiniVolcano() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const [pts, setPts] = useState<VolcanoPoint[] | null>(null);
  const [failed, setFailed] = useState(false);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("volcano_preview.json"))
      .then((r) => (r.ok ? r.json() : Promise.reject(r.status)))
      .then((d: VolcanoPoint[]) => {
        if (!alive) return;
        // Draw non-significant first so significant marks sit on top.
        d.sort((a, b) => Number(a.sig) - Number(b.sig));
        setPts(d);
      })
      .catch(() => alive && setFailed(true));
    return () => {
      alive = false;
    };
  }, []);

  if (failed) return null;

  const ready = pts && pts.length > 0 && width > 0 && height > 0;

  return (
    <div ref={ref} className="h-full w-full">
      {!ready ? (
        <div className="h-full w-full animate-pulse rounded-md bg-muted/40 motion-reduce:animate-none" />
      ) : (
        (() => {
          const xMax = Math.max(1e-3, ...pts!.map((p) => Math.abs(p.x)));
          const yMax = Math.max(1e-3, ...pts!.map((p) => p.y));
          const innerW = width - 2 * PAD;
          const innerH = height - 2 * PAD;
          const sx = (x: number) => PAD + ((x + xMax) / (2 * xMax)) * innerW;
          const sy = (y: number) => height - PAD - (y / yMax) * innerH;
          const r = Math.max(1.1, Math.min(2.4, width / 130));
          const zeroX = sx(0);
          return (
            <svg
              width={width}
              height={height}
              className="block"
              role="img"
              aria-label="Volcano plot of differential expression"
            >
              <line
                x1={zeroX}
                x2={zeroX}
                y1={PAD}
                y2={height - PAD}
                stroke="var(--color-border)"
                strokeWidth={1}
                strokeDasharray="2 3"
                shapeRendering="crispEdges"
              />
              {pts!.map((p, i) => {
                const fill = p.sig
                  ? p.x > 0
                    ? DIVERGING.up
                    : DIVERGING.down
                  : CONTROL;
                return (
                  <circle
                    key={i}
                    cx={sx(p.x)}
                    cy={sy(p.y)}
                    r={r}
                    fill={fill}
                    fillOpacity={p.sig ? SIG_OPACITY : NONSIG_OPACITY}
                  />
                );
              })}
            </svg>
          );
        })()
      )}
    </div>
  );
}
