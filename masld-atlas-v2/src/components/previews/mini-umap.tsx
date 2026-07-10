"use client";

/**
 * Landing mini-viz: a thumbnail UMAP of the single-cell atlas.
 *
 * Self-contained and fail-soft — fetches `umap_thumbnail.json`
 * (`[{x,y,celltype}]`, ~2.5k points) and paints them on a
 * devicePixelRatio-aware canvas, coloring each point by cell type via the
 * colorblind-safe categorical palette (a stable celltype -> color map). Canvas
 * marks read palette hex directly (canvas can't resolve theme CSS vars); the
 * saturated categorical hues read on both light and dark card surfaces. A tiny
 * DOM legend labels the most abundant cell types.
 */

import { useEffect, useMemo, useState } from "react";
import { dataUrl } from "@/lib/data-base";
import { cellTypeColor } from "@/lib/palette";
import { useMeasure } from "./use-measure";

interface UmapPoint {
  x: number;
  y: number;
  celltype: string;
}

const PAD = 4;
const DOT = 1.5;
const DOT_ALPHA = 0.6;

export function MiniUmap() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const [pts, setPts] = useState<UmapPoint[] | null>(null);
  const [failed, setFailed] = useState(false);
  const [canvas, setCanvas] = useState<HTMLCanvasElement | null>(null);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("umap_thumbnail.json"))
      .then((r) => (r.ok ? r.json() : Promise.reject(r.status)))
      .then((d: UmapPoint[]) => alive && setPts(d))
      .catch(() => alive && setFailed(true));
    return () => {
      alive = false;
    };
  }, []);

  // Stable celltype -> color map (alphabetical order → deterministic).
  const colorMap = useMemo(() => {
    const m = new Map<string, string>();
    if (!pts) return m;
    [...new Set(pts.map((p) => p.celltype))]
      .sort()
      .forEach((ct, i) => m.set(ct, cellTypeColor(ct, i)));
    return m;
  }, [pts]);

  const legend = useMemo(() => {
    if (!pts) return [];
    const counts = new Map<string, number>();
    for (const p of pts) counts.set(p.celltype, (counts.get(p.celltype) ?? 0) + 1);
    return [...counts.entries()]
      .sort((a, b) => b[1] - a[1])
      .slice(0, 4)
      .map(([ct]) => ct);
  }, [pts]);

  useEffect(() => {
    if (!canvas || !pts || pts.length === 0 || width === 0 || height === 0) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);

    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of pts) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.y > maxY) maxY = p.y;
    }
    const spanX = maxX - minX || 1;
    const spanY = maxY - minY || 1;
    const iw = width - 2 * PAD;
    const ih = height - 2 * PAD;
    ctx.globalAlpha = DOT_ALPHA;
    for (const p of pts) {
      const px = PAD + ((p.x - minX) / spanX) * iw;
      // Flip Y so higher UMAP-2 is toward the top.
      const py = PAD + (1 - (p.y - minY) / spanY) * ih;
      ctx.fillStyle = colorMap.get(p.celltype) ?? "#9E9E9E";
      ctx.fillRect(px - DOT / 2, py - DOT / 2, DOT, DOT);
    }
    ctx.globalAlpha = 1;
  }, [canvas, pts, colorMap, width, height]);

  if (failed) return null;

  const ready = pts && pts.length > 0;

  return (
    <div ref={ref} className="relative h-full w-full">
      {!ready ? (
        <div className="h-full w-full animate-pulse rounded-md bg-muted/40 motion-reduce:animate-none" />
      ) : (
        <>
          <canvas ref={setCanvas} className="block" aria-label="Single-cell UMAP thumbnail" />
          {legend.length > 0 && (
            <ul className="pointer-events-none absolute left-1 top-1 flex flex-col gap-0.5">
              {legend.map((ct) => (
                <li
                  key={ct}
                  className="flex items-center gap-1 text-[9px] leading-none text-muted-foreground"
                >
                  <span
                    aria-hidden
                    className="size-1.5 shrink-0 rounded-full"
                    style={{ backgroundColor: colorMap.get(ct) }}
                  />
                  <span className="max-w-[90px] truncate">{ct}</span>
                </li>
              ))}
            </ul>
          )}
        </>
      )}
    </div>
  );
}
