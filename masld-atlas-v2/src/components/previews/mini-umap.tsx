"use client";

/**
 * Landing mini-viz: an animated, glowing UMAP of the single-cell atlas.
 *
 * Fetches `umap_thumbnail.json` (`[{x,y,celltype}]`, ~2.5k points) and paints
 * them on a devicePixelRatio-aware canvas as soft ROUND dots, colored by cell
 * type via the canonical Fig-3 palette (`cellTypeColor`). On dark (the default
 * theme) dots are drawn with additive blending so dense clusters bloom like a
 * nebula; on light they are soft solid dots (additive would blow out to white).
 *
 * Motion: points fade + scale into place (a ~1.2s "bloom-in"), then breathe
 * with a very slow ambient drift. Perf: one cached glow sprite per cell-type
 * color is blitted per point (fast), the drift is throttled to ~30fps, the
 * loop pauses when the tile scrolls offscreen (IntersectionObserver), and
 * reduced motion collapses to a single static frame. A tiny DOM legend labels
 * the most abundant cell types.
 */

import { useEffect, useMemo, useRef, useState } from "react";
import { dataUrl } from "@/lib/data-base";
import { cellTypeColor } from "@/lib/palette";
import { useMeasure } from "./use-measure";

interface UmapPoint {
  x: number;
  y: number;
  celltype: string;
}

const PAD = 6;
const BLOOM_MS = 1200;
const STAGGER_MS = 340; // spread of per-point bloom start
const DRIFT_AMP = 1.5; // px
const DRIFT_FPS = 30;

/** "#RRGGBB" -> "rgba(r,g,b,a)". */
function rgba(hex: string, a: number): string {
  const h = hex.replace("#", "");
  const n = parseInt(
    h.length === 3 ? h.replace(/(.)/g, "$1$1") : h.slice(0, 6),
    16
  );
  return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${a})`;
}

export function MiniUmap() {
  const [ref, { width, height }] = useMeasure<HTMLDivElement>();
  const [pts, setPts] = useState<UmapPoint[] | null>(null);
  const [failed, setFailed] = useState(false);
  const [canvas, setCanvas] = useState<HTMLCanvasElement | null>(null);
  const [dark, setDark] = useState(true);
  const visibleRef = useRef(true);

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

  // Track the active theme so the glow (additive on dark) recolors on toggle.
  useEffect(() => {
    const read = () =>
      setDark(document.documentElement.classList.contains("dark"));
    read();
    const obs = new MutationObserver(read);
    obs.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });
    return () => obs.disconnect();
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

    const reduce = window.matchMedia?.(
      "(prefers-reduced-motion: reduce)"
    ).matches;

    // Dot radius scales with the tile so the 2x2 hero reads bigger than a 1x1.
    const minDim = Math.min(width, height);
    const r = Math.max(1.8, Math.min(3.2, minDim / 110));
    const spriteR = Math.ceil(r * (dark ? 3.4 : 2.0)); // halo radius (css px)

    // Screen positions, per-point phase seeds, and colors.
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
    const N = pts.length;
    const sx = new Float32Array(N);
    const sy = new Float32Array(N);
    const seed = new Float32Array(N);
    const colors: string[] = new Array(N);
    for (let i = 0; i < N; i++) {
      const p = pts[i];
      sx[i] = PAD + ((p.x - minX) / spanX) * iw;
      // Flip Y so higher UMAP-2 is toward the top.
      sy[i] = PAD + (1 - (p.y - minY) / spanY) * ih;
      seed[i] = (i * 0.6180339887) % 1; // deterministic 0..1
      colors[i] = colorMap.get(p.celltype) ?? "#9E9E9E";
    }

    // One cached glow sprite per color; blitting is far cheaper than per-point
    // radial gradients every frame.
    const sprites = new Map<string, HTMLCanvasElement>();
    const spriteFor = (color: string): HTMLCanvasElement => {
      const hit = sprites.get(color);
      if (hit) return hit;
      const s = document.createElement("canvas");
      s.width = s.height = Math.ceil(spriteR * 2 * dpr);
      const sc = s.getContext("2d")!;
      sc.scale(dpr, dpr);
      const g = sc.createRadialGradient(spriteR, spriteR, 0, spriteR, spriteR, spriteR);
      if (dark) {
        g.addColorStop(0, rgba(color, 0.95));
        g.addColorStop(0.4, rgba(color, 0.5));
        g.addColorStop(1, rgba(color, 0));
      } else {
        g.addColorStop(0, rgba(color, 0.9));
        g.addColorStop(0.55, rgba(color, 0.55));
        g.addColorStop(1, rgba(color, 0));
      }
      sc.fillStyle = g;
      sc.beginPath();
      sc.arc(spriteR, spriteR, spriteR, 0, Math.PI * 2);
      sc.fill();
      sprites.set(color, s);
      return s;
    };
    for (const c of new Set(colors)) spriteFor(c);

    const composite: GlobalCompositeOperation = dark ? "lighter" : "source-over";

    const paint = (bloom: number, now: number, drift: boolean) => {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, width, height);
      ctx.globalCompositeOperation = composite;
      for (let i = 0; i < N; i++) {
        // Per-point staggered bloom (fade + scale in).
        const delay = seed[i] * STAGGER_MS;
        const lt = Math.max(0, Math.min(1, (bloom - delay) / (BLOOM_MS - STAGGER_MS)));
        const ease = 1 - Math.pow(1 - lt, 3);
        if (ease <= 0) continue;
        let ox = 0, oy = 0;
        if (drift) {
          const ph = seed[i] * Math.PI * 2;
          ox = Math.cos(now * 0.00016 + ph) * DRIFT_AMP;
          oy = Math.sin(now * 0.00014 + ph) * DRIFT_AMP;
        }
        const drawR = spriteR * (0.28 + 0.72 * ease);
        ctx.globalAlpha = ease;
        ctx.drawImage(spriteFor(colors[i]), sx[i] + ox - drawR, sy[i] + oy - drawR, drawR * 2, drawR * 2);
      }
      ctx.globalAlpha = 1;
      ctx.globalCompositeOperation = "source-over";
    };

    if (reduce) {
      paint(BLOOM_MS, 0, false); // static final frame, no motion
      return;
    }

    let raf = 0;
    let start = 0;
    let lastDrift = 0;
    const driftInterval = 1000 / DRIFT_FPS;
    const frame = (now: number) => {
      if (!start) start = now;
      const bloom = now - start;
      const blooming = bloom < BLOOM_MS;
      if (blooming) {
        paint(bloom, now, false); // full fps while blooming
      } else if (now - lastDrift >= driftInterval) {
        lastDrift = now;
        paint(BLOOM_MS, now, true); // throttled ambient drift
      }
      if (blooming || visibleRef.current) raf = requestAnimationFrame(frame);
    };
    raf = requestAnimationFrame(frame);

    // Pause the loop while the tile is offscreen; resume (drift only) on return.
    const io = new IntersectionObserver(
      (entries) => {
        const vis = entries[0]?.isIntersecting ?? true;
        const was = visibleRef.current;
        visibleRef.current = vis;
        if (vis && !was) {
          cancelAnimationFrame(raf);
          raf = requestAnimationFrame(frame);
        }
      },
      { threshold: 0 }
    );
    io.observe(canvas);

    return () => {
      cancelAnimationFrame(raf);
      io.disconnect();
    };
  }, [canvas, pts, colorMap, width, height, dark]);

  if (failed) return null;

  const ready = pts && pts.length > 0;

  return (
    <div ref={ref} className="relative h-full w-full">
      {!ready ? (
        <div className="h-full w-full animate-pulse rounded-md bg-muted/40 motion-reduce:animate-none" />
      ) : (
        <>
          <canvas
            ref={setCanvas}
            className="block"
            aria-label="Single-cell UMAP thumbnail"
            style={{
              maskImage:
                "radial-gradient(135% 135% at 50% 45%, black 58%, transparent 100%)",
              WebkitMaskImage:
                "radial-gradient(135% 135% at 50% 45%, black 58%, transparent 100%)",
            }}
          />
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
