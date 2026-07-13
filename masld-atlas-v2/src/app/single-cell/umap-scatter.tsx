"use client";

/**
 * WebGL UMAP scatter over the downsampled single-cell embedding
 * (`sc_umap_downsampled.parquet`, ~60k cells) via regl-scatterplot.
 *
 * This module is browser-only — the page loads it through
 * `dynamic(() => import("./umap-scatter"), { ssr: false })`, so the top-level
 * regl-scatterplot import never runs on the server.
 *
 * Coloring is by cell type via the canonical `cellTypeColor` map (matching the
 * published Fig 3 palette) OR by disease state (Healthy → control gray, disease →
 * diverging red). A "reduce points" toggle strides the drawn set for lower-power
 * devices; a "lasso" toggle flips the drag behaviour to marquee-select, whose
 * `select` / `deselect` events lift the selected cells to the page for a live
 * per-cell-type / condition summary. Recolor/reduce redraws tween positions via
 * `draw({ transition })`, gated on reduced motion.
 */

import { useEffect, useMemo, useRef, useState } from "react";
import createScatterplot from "regl-scatterplot";
import { useReducedMotion } from "framer-motion";
import { CONTROL, DIVERGING, cellTypeColor } from "@/lib/palette";

export interface UmapPoint {
  x: number;
  y: number;
  cellType: string;
  dataset: string;
  condition: string;
  cellId: string;
}

export type ColorMode = "cell_type" | "condition_harmonized";

interface Props {
  points: UmapPoint[];
  colorMode: ColorMode;
  reduce: boolean;
  /** When true, a plain drag draws a marquee lasso instead of panning. */
  lasso: boolean;
  /** Selected cells (via lasso) lifted to the page; [] on deselect. */
  onSelect?: (pts: UmapPoint[]) => void;
  /** Bump to imperatively clear the on-canvas lasso selection. */
  clearNonce?: number;
  height?: number;
}

/** Build the category index + ordered color palette for the active mode. */
function encoding(points: UmapPoint[], mode: ColorMode) {
  if (mode === "condition_harmonized") {
    // Healthy first (→ control gray), then everything else (disease red).
    const cats = Array.from(new Set(points.map((p) => p.condition))).sort((a, b) =>
      a === "Healthy" ? -1 : b === "Healthy" ? 1 : a.localeCompare(b)
    );
    const colors = cats.map((c) => (/^healthy$/i.test(c) ? CONTROL : DIVERGING.up));
    const index = new Map(cats.map((c, i) => [c, i]));
    return { cats, colors, keyOf: (p: UmapPoint) => index.get(p.condition) ?? 0 };
  }
  const cats = Array.from(new Set(points.map((p) => p.cellType))).sort();
  const colors = cats.map((c, i) => cellTypeColor(c, i));
  const index = new Map(cats.map((c, i) => [c, i]));
  return { cats, colors, keyOf: (p: UmapPoint) => index.get(p.cellType) ?? 0 };
}

export default function UmapScatter({
  points,
  colorMode,
  reduce,
  lasso,
  onSelect,
  clearNonce = 0,
  height = 520,
}: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const spRef = useRef<ReturnType<typeof createScatterplot> | null>(null);
  // Metadata aligned to the CURRENTLY drawn points (order matches the draw call).
  const drawnMetaRef = useRef<UmapPoint[]>([]);
  // Keep the select callback in a ref so the one-time subscription stays fresh.
  const onSelectRef = useRef<Props["onSelect"]>(onSelect);
  useEffect(() => {
    onSelectRef.current = onSelect;
  });
  const reduceMotion = useReducedMotion();
  const [hover, setHover] = useState<UmapPoint | null>(null);

  // Normalize UMAP coords to [-1, 1] once, preserving aspect ratio.
  const normalized = useMemo(() => {
    if (points.length === 0) return [] as { nx: number; ny: number }[];
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of points) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.y > maxY) maxY = p.y;
    }
    const cx = (minX + maxX) / 2;
    const cy = (minY + maxY) / 2;
    const span = Math.max(maxX - minX, maxY - minY) / 2 || 1;
    return points.map((p) => ({ nx: (p.x - cx) / span, ny: (p.y - cy) / span }));
  }, [points]);

  const { cats, colors, keyOf } = useMemo(
    () => encoding(points, colorMode),
    [points, colorMode]
  );

  // Create the scatterplot once.
  useEffect(() => {
    const canvas = canvasRef.current;
    const container = containerRef.current;
    if (!canvas || !container) return;
    const sp = createScatterplot({
      canvas,
      width: container.clientWidth || 640,
      height,
      pointSize: 2,
      opacity: 0.68,
      backgroundColor: [0, 0, 0, 0],
      lassoColor: [0.62, 0.62, 0.62, 1], // control gray, chrome only
    });
    spRef.current = sp;

    const over = sp.subscribe("pointOver", (i: number) => {
      const m = drawnMetaRef.current[i];
      if (m) setHover(m);
    });
    const out = sp.subscribe("pointOut", () => setHover(null));
    const sel = sp.subscribe("select", ({ points: idxs }: { points: number[] }) => {
      const picked = idxs
        .map((i) => drawnMetaRef.current[i])
        .filter((m): m is UmapPoint => Boolean(m));
      onSelectRef.current?.(picked);
    });
    const desel = sp.subscribe("deselect", () => onSelectRef.current?.([]));

    const onResize = () => {
      const w = container.clientWidth;
      if (w > 0) sp.set({ width: w, height });
    };
    window.addEventListener("resize", onResize);

    return () => {
      window.removeEventListener("resize", onResize);
      sp.unsubscribe(over);
      sp.unsubscribe(out);
      sp.unsubscribe(sel);
      sp.unsubscribe(desel);
      sp.destroy();
      spRef.current = null;
    };
  }, [height]);

  // Toggle marquee-lasso vs pan/zoom on the drag gesture.
  useEffect(() => {
    const sp = spRef.current;
    if (!sp) return;
    sp.set({ mouseMode: lasso ? "lasso" : "panZoom" });
  }, [lasso]);

  // Imperative clear of the on-canvas selection (page "Clear" button). Skips
  // the initial mount so it doesn't fire a spurious deselect.
  const firstClear = useRef(true);
  useEffect(() => {
    if (firstClear.current) {
      firstClear.current = false;
      return;
    }
    spRef.current?.deselect();
  }, [clearNonce]);

  // (Re)draw whenever the data, coloring, or reduce toggle changes.
  useEffect(() => {
    const sp = spRef.current;
    if (!sp || normalized.length === 0) return;
    const stride = reduce ? 2 : 1;
    const drawn: number[][] = [];
    const meta: UmapPoint[] = [];
    for (let i = 0; i < points.length; i += stride) {
      drawn.push([normalized[i].nx, normalized[i].ny, keyOf(points[i])]);
      meta.push(points[i]);
    }
    drawnMetaRef.current = meta;
    sp.set({ pointColor: colors, colorBy: "category", pointSize: reduce ? 3 : 2 });
    void sp.draw(drawn, { transition: !reduceMotion, transitionDuration: 500 });
  }, [points, normalized, colors, keyOf, reduce, reduceMotion]);

  return (
    <div className="space-y-3">
      <div
        ref={containerRef}
        className="relative w-full overflow-hidden rounded-lg border border-border bg-card"
        style={{ height }}
      >
        <canvas ref={canvasRef} className="block h-full w-full" />
        {hover && (
          <div className="pointer-events-none absolute right-3 top-3 rounded-md border border-border bg-background/90 px-3 py-2 text-xs shadow-sm backdrop-blur">
            <div className="font-semibold">{hover.cellType}</div>
            <div className="text-muted-foreground">{hover.condition}</div>
            <div className="font-mono text-[10px] text-muted-foreground">
              {hover.dataset}
            </div>
          </div>
        )}
        {lasso && (
          <div className="pointer-events-none absolute left-3 top-3 rounded-md border border-border bg-background/90 px-2 py-1 text-[10px] text-muted-foreground shadow-sm backdrop-blur">
            Drag to lasso · Esc / click empty to clear
          </div>
        )}
      </div>
      {/* Legend */}
      <ul className="flex flex-wrap gap-x-4 gap-y-1.5 text-xs text-muted-foreground">
        {cats.map((c, i) => (
          <li key={c} className="flex items-center gap-1.5">
            <span
              aria-hidden
              className="size-2.5 shrink-0 rounded-[3px]"
              style={{ backgroundColor: colors[i] }}
            />
            <span>{c}</span>
          </li>
        ))}
      </ul>
    </div>
  );
}
