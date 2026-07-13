"use client";

/**
 * Live gene-constellation for the landing hero — DECORATIVE CHROME, not a data
 * figure. Renders the top-genes subgraph from `hero-data.ts` on a transparent
 * <canvas> via react-force-graph-2d: brand-gradient "stars" sized by evidence,
 * faint pathway/drug connector hubs, faint links. No interaction, no labels, no
 * quantitative color scale — it settles for a fixed number of ticks then
 * freezes. Colors are resolved from the brand CSS tokens (`--primary`,
 * `--accent-brand`, `--muted-foreground`) so it tracks the light/dark theme.
 *
 * Dynamically imported with `ssr:false` (needs a browser canvas). The hero only
 * mounts this on an idle callback when motion is allowed and canvas is capable,
 * so the static SVG fallback covers first paint / reduced-motion / no-WebGL.
 */

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { HeroGraph, HeroNode } from "@/lib/hero-data";

// Dynamic import to bypass SSR — react-force-graph requires <canvas>.
const ForceGraph2D = dynamic(() => import("react-force-graph-2d"), {
  ssr: false,
});

const TWO_PI = 2 * Math.PI;

// ---------------------------------------------------------------------------
// Token color resolution (decorative chrome colors — resolved to rgb so the
// canvas can paint them and interpolate; theme-reactive via a class observer).
// ---------------------------------------------------------------------------

type Rgb = [number, number, number];

interface Colors {
  primary: Rgb;
  accent: Rgb;
  hub: Rgb;
}

const FALLBACK: Colors = {
  primary: [88, 136, 252],
  accent: [150, 96, 230],
  hub: [130, 130, 140],
};

function parseToRgb(color: string): Rgb | null {
  if (typeof document === "undefined") return null;
  const cv = document.createElement("canvas");
  cv.width = cv.height = 1;
  const ctx = cv.getContext("2d");
  if (!ctx) return null;
  try {
    ctx.fillStyle = "#000";
    ctx.fillStyle = color; // invalid strings leave fillStyle at the prior value
    ctx.fillRect(0, 0, 1, 1);
    const d = ctx.getImageData(0, 0, 1, 1).data;
    return [d[0], d[1], d[2]];
  } catch {
    return null;
  }
}

function readVar(name: string): Rgb | null {
  if (typeof document === "undefined") return null;
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return v ? parseToRgb(v) : null;
}

function resolveColors(): Colors {
  if (typeof document === "undefined") return FALLBACK;
  return {
    primary: readVar("--primary") ?? FALLBACK.primary,
    accent: readVar("--accent-brand") ?? FALLBACK.accent,
    hub: readVar("--muted-foreground") ?? FALLBACK.hub,
  };
}

const lerp = (a: Rgb, b: Rgb, t: number): Rgb => [
  a[0] + (b[0] - a[0]) * t,
  a[1] + (b[1] - a[1]) * t,
  a[2] + (b[2] - a[2]) * t,
];

const css = (c: Rgb, a = 1): string =>
  `rgba(${Math.round(c[0])}, ${Math.round(c[1])}, ${Math.round(c[2])}, ${a})`;

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

interface PaintNode extends HeroNode {
  x?: number;
  y?: number;
  __r: number;
}

export function GeneConstellation({ graph }: { graph: HeroGraph }) {
  const containerRef = useRef<HTMLDivElement>(null);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const fgRef = useRef<any>(null);
  const colorsRef = useRef<Colors>(resolveColors());
  const [dim, setDim] = useState({ w: 0, h: 0 });

  // Clone into fresh objects — react-force-graph mutates node/link entries
  // (adds x/y, replaces source/target with node refs); never mutate the cache.
  const graphData = useMemo(() => {
    const nodes = graph.nodes.map((n) => ({
      ...n,
      __r: n.kind === "gene" ? 3.2 + n.weight * 6 : 2.2,
    }));
    const links = graph.links.map((l) => ({ source: l.source, target: l.target }));
    return { nodes, links };
  }, [graph]);

  // Measure container.
  useEffect(() => {
    const el = containerRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      for (const e of entries) {
        const { width, height } = e.contentRect;
        if (width > 0 && height > 0)
          setDim({ w: Math.floor(width), h: Math.floor(height) });
      }
    });
    ro.observe(el);
    const rect = el.getBoundingClientRect();
    if (rect.width > 0 && rect.height > 0)
      setDim({ w: Math.floor(rect.width), h: Math.floor(rect.height) });
    return () => ro.disconnect();
  }, []);

  // Re-resolve colors on theme change and repaint the (frozen) canvas.
  useEffect(() => {
    if (typeof document === "undefined") return;
    const update = () => {
      colorsRef.current = resolveColors();
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (fgRef.current as any)?.d3ReheatSimulation?.();
    };
    update();
    const obs = new MutationObserver(update);
    obs.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });
    return () => obs.disconnect();
  }, []);

  // Loosen forces for an airy, spread-out constellation.
  useEffect(() => {
    const fg = fgRef.current;
    if (!fg) return;
    fg.d3Force("charge")?.strength(-70);
    fg.d3Force("link")?.distance(40);
  }, [graphData]);

  const paintNode = useCallback(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (node: any, ctx: CanvasRenderingContext2D) => {
      const n = node as PaintNode;
      if (n.x == null || n.y == null) return;
      const c = colorsRef.current;
      const isGene = n.kind === "gene";
      const fill = isGene ? lerp(c.primary, c.accent, n.weight) : c.hub;
      ctx.save();
      ctx.globalAlpha = isGene ? 0.95 : 0.5;
      if (isGene) {
        ctx.shadowBlur = 14;
        ctx.shadowColor = css(fill, 0.9);
      }
      ctx.beginPath();
      ctx.arc(n.x, n.y, n.__r, 0, TWO_PI);
      ctx.fillStyle = css(fill);
      ctx.fill();
      ctx.restore();
    },
    []
  );

  const linkColor = useCallback(() => css(colorsRef.current.primary, 0.26), []);

  const onEngineStop = useCallback(() => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (fgRef.current as any)?.zoomToFit?.(500, 26);
  }, []);

  return (
    <div ref={containerRef} className="h-full w-full">
      {dim.w > 0 && dim.h > 0 && (
        <ForceGraph2D
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
          ref={fgRef as any}
          graphData={graphData}
          width={dim.w}
          height={dim.h}
          nodeCanvasObject={paintNode}
          linkColor={linkColor}
          linkWidth={1}
          cooldownTicks={70}
          d3AlphaDecay={0.03}
          enableNodeDrag={false}
          enableZoomInteraction={false}
          enablePanInteraction={false}
          enablePointerInteraction={false}
          backgroundColor="transparent"
          onEngineStop={onEngineStop}
        />
      )}
    </div>
  );
}

export default GeneConstellation;
