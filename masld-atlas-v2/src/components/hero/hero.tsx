"use client";

/**
 * Immersive landing hero — decorative brand chrome. A full-bleed brand-gradient
 * + vignette backdrop with a live gene-constellation behind a display-scale
 * gradient headline and the primary CTAs.
 *
 * The static SVG fallback paints first (correct with JS disabled / no canvas /
 * reduced motion / Save-Data); the live canvas hydrates on an idle callback
 * only when motion is allowed and the device is capable. Nothing here encodes
 * data — no logFC scale, no legend.
 */

import { useEffect, useState } from "react";
import { HashLink as Link } from "@/components/hash-link";
import dynamic from "next/dynamic";
import { Button } from "@/components/ui/button";
import { useAppStore } from "@/lib/store";
import { loadHeroGraph, type HeroGraph } from "@/lib/hero-data";
import { HeroFallback } from "@/components/hero/hero-fallback";
import { Aurora } from "@/components/backgrounds/aurora";

const GeneConstellation = dynamic(
  () => import("@/components/hero/gene-constellation"),
  { ssr: false, loading: () => <HeroFallback /> }
);

export function Hero() {
  const { setCommandOpen } = useAppStore();
  const [graph, setGraph] = useState<HeroGraph | null>(null);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const reduce = window.matchMedia?.(
      "(prefers-reduced-motion: reduce)"
    ).matches;
    const canvasOk = !!document.createElement("canvas").getContext?.("2d");
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const saveData = (navigator as any).connection?.saveData;
    if (reduce || !canvasOk || saveData) return; // keep the static fallback

    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const w = window as any;
    let handle: number;
    const run = () => {
      loadHeroGraph(40)
        .then((g) => {
          if (g.nodes.length > 0) setGraph(g);
        })
        .catch(() => {});
    };
    if (typeof w.requestIdleCallback === "function") {
      handle = w.requestIdleCallback(run, { timeout: 1800 });
    } else {
      handle = window.setTimeout(run, 400);
    }
    return () => {
      if (typeof w.cancelIdleCallback === "function") w.cancelIdleCallback(handle);
      else window.clearTimeout(handle);
    };
  }, []);

  return (
    <section className="relative isolate overflow-hidden">
      {/* Backmost ambient aurora — decorative brand glow beneath every layer. */}
      <Aurora />
      {/* Live constellation (or static fallback) behind everything. */}
      <div
        aria-hidden
        className="pointer-events-none absolute inset-0 -z-10 opacity-70 dark:opacity-80"
        style={{
          maskImage:
            "radial-gradient(115% 92% at 50% 42%, transparent 0%, transparent 30%, black 70%)",
          WebkitMaskImage:
            "radial-gradient(115% 92% at 50% 42%, transparent 0%, transparent 30%, black 70%)",
        }}
      >
        {graph ? <GeneConstellation graph={graph} /> : <HeroFallback />}
      </div>
      {/* Brand glow + depth vignette layered over the constellation. */}
      <div
        aria-hidden
        className="pointer-events-none absolute inset-0 -z-10"
        style={{ backgroundImage: "var(--gradient-hero)" }}
      />
      <div
        aria-hidden
        className="pointer-events-none absolute inset-0 -z-10"
        style={{ backgroundImage: "var(--depth-vignette)" }}
      />

      <div className="mx-auto max-w-3xl px-6 pt-14 pb-16 text-center sm:pt-20 sm:pb-20">
        <h1 className="text-gradient font-display text-display-xl font-semibold tracking-tight text-balance">
          MASLD Atlas
        </h1>
        <p className="mx-auto mt-5 max-w-xl text-lg text-muted-foreground">
          A multi-modal atlas for metabolic dysfunction-associated steatotic
          liver disease — genetics, transcriptomics, and drug evidence in one
          place.
        </p>
        <div className="mt-8 flex items-center justify-center gap-3">
          <Button size="lg" onClick={() => setCommandOpen(true)}>
            <svg
              className="size-4"
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
              strokeWidth={2}
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M21 21l-5.197-5.197m0 0A7.5 7.5 0 105.196 5.196a7.5 7.5 0 0010.607 10.607z"
              />
            </svg>
            Search genes
          </Button>
          <Button variant="outline" size="lg" render={<Link href="/atlas" />}>
            Explore atlas
          </Button>
        </div>
      </div>
    </section>
  );
}

export default Hero;
