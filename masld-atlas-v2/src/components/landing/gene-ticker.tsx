"use client";

/**
 * <GeneTicker> — a slow marquee of the top differentially expressed genes.
 *
 * Fetches `landing_ticker.json` (`{symbol, logFC, dir}[]`) and renders each as a
 * small pill chip linking to its gene page. The up/down arrow COLOR is the only
 * colored element — a legitimate data encoding of regulation direction from the
 * DIVERGING scale (up = red, down = blue). Gene text stays uncolored.
 *
 * Fail-soft: renders nothing if the fetch fails or returns an empty list.
 * Offscreen-pause: an IntersectionObserver freezes the scroll when the row
 * leaves the viewport (passed to <MarqueeRow paused>).
 */

import { useEffect, useRef, useState } from "react";
import { ArrowDown, ArrowUp } from "lucide-react";
import { HashLink } from "@/components/hash-link";
import { MarqueeRow } from "@/components/ui/marquee";
import { dataUrl } from "@/lib/data-base";
import { DIVERGING } from "@/lib/palette";
import { cn } from "@/lib/utils";

interface TickerGene {
  symbol: string;
  logFC: number;
  dir: "up" | "down";
}

export function GeneTicker({ className }: { className?: string }) {
  const [genes, setGenes] = useState<TickerGene[]>([]);
  const [paused, setPaused] = useState(false);
  const rootRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let alive = true;
    fetch(dataUrl("landing_ticker.json"))
      .then((r) => r.json())
      .then((data: TickerGene[]) => {
        if (alive && Array.isArray(data)) setGenes(data);
      })
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, []);

  // Freeze the scroll while the ticker is offscreen (saves paint work).
  useEffect(() => {
    const el = rootRef.current;
    if (!el || typeof IntersectionObserver === "undefined") return;
    const io = new IntersectionObserver(
      ([entry]) => setPaused(!entry.isIntersecting),
      { rootMargin: "0px" }
    );
    io.observe(el);
    return () => io.disconnect();
  }, [genes.length]);

  // Fail-soft: nothing to show (fetch failed / empty / not yet loaded).
  if (genes.length === 0) return null;

  return (
    <div
      ref={rootRef}
      aria-label="Top differentially expressed genes"
      className={className}
    >
      <MarqueeRow durationSec={60} paused={paused}>
        {genes.map((g, i) => {
          const Arrow = g.dir === "up" ? ArrowUp : ArrowDown;
          const color = g.dir === "up" ? DIVERGING.up : DIVERGING.down;
          return (
            <HashLink
              key={`${g.symbol}-${i}`}
              href={`#/gene?symbol=${g.symbol}`}
              className={cn(
                "mx-1.5 inline-flex items-center gap-1 rounded-full border border-border",
                "bg-card px-2.5 py-1 text-xs transition-colors hover:border-primary/40"
              )}
            >
              <span className="font-mono font-medium tabular-nums">
                {g.symbol}
              </span>
              <Arrow className="size-3 shrink-0" style={{ color }} aria-hidden />
              <span className="sr-only">
                {g.dir === "up" ? "upregulated" : "downregulated"}
              </span>
            </HashLink>
          );
        })}
      </MarqueeRow>
    </div>
  );
}

export default GeneTicker;
