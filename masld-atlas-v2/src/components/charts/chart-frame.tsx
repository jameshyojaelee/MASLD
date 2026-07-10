"use client";

import { useEffect, useRef, useState, type ReactNode } from "react";
import { Group } from "@visx/group";
import { cn } from "@/lib/utils";

export interface ChartMargin {
  top: number;
  right: number;
  bottom: number;
  left: number;
}

const DEFAULT_MARGIN: ChartMargin = { top: 8, right: 12, bottom: 32, left: 44 };

export interface ChartFrameRenderArgs {
  /** Plotting area width (svg width minus horizontal margins). */
  innerWidth: number;
  /** Plotting area height (svg height minus vertical margins). */
  innerHeight: number;
  margin: ChartMargin;
}

export interface ChartFrameProps {
  /** Optional heading rendered above the plot (theme foreground, never colored). */
  title?: ReactNode;
  /** Optional caption rendered below the plot in muted text. */
  caption?: ReactNode;
  /** Fixed height of the SVG in px. Width is measured from the container. */
  height?: number;
  /** Explicit width override; when omitted the container is measured. */
  width?: number;
  margin?: Partial<ChartMargin>;
  className?: string;
  /** Accessible label for the SVG. */
  ariaLabel?: string;
  /**
   * Render-prop receiving the resolved plotting-area dimensions. Content is
   * already translated by the margins (wrapped in a visx <Group>), so draw in
   * [0, innerWidth] × [0, innerHeight].
   */
  children: (args: ChartFrameRenderArgs) => ReactNode;
}

/**
 * Responsive SVG chart shell shared by every visx chart in the portal.
 *
 * Handles the boilerplate — container measurement, margin convention, the
 * translated <Group>, and consistent title/caption chrome — so individual
 * charts only implement their marks. This is intentionally a thin skeleton;
 * concrete charts (Volcano, ForestPlot, Heatmap, …) build on top of it.
 */
export function ChartFrame({
  title,
  caption,
  height = 320,
  width: widthProp,
  margin: marginProp,
  className,
  ariaLabel,
  children,
}: ChartFrameProps) {
  const margin = { ...DEFAULT_MARGIN, ...marginProp };
  const containerRef = useRef<HTMLDivElement>(null);
  const [measured, setMeasured] = useState(0);

  useEffect(() => {
    if (widthProp != null) return;
    const el = containerRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width ?? 0;
      setMeasured(w);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [widthProp]);

  const width = widthProp ?? measured;
  const innerWidth = Math.max(0, width - margin.left - margin.right);
  const innerHeight = Math.max(0, height - margin.top - margin.bottom);
  const ready = width > 0;

  // First-paint fade: once the container is measured (width > 0) the marks
  // cross-fade in instead of popping. Skipped entirely under reduced motion.
  const [painted, setPainted] = useState(false);
  const reduceMotion =
    typeof window !== "undefined" &&
    typeof window.matchMedia === "function" &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  useEffect(() => {
    if (ready) setPainted(true);
  }, [ready]);

  return (
    <figure className={cn("m-0 w-full", className)}>
      {title && (
        <figcaption className="mb-1 text-sm font-medium text-foreground">
          {title}
        </figcaption>
      )}
      <div ref={containerRef} className="w-full">
        {ready && (
          <svg
            width={width}
            height={height}
            role="img"
            aria-label={ariaLabel ?? (typeof title === "string" ? title : undefined)}
            style={
              reduceMotion
                ? undefined
                : {
                    opacity: painted ? 1 : 0,
                    transition: "opacity var(--duration-base) var(--ease-out-expo)",
                  }
            }
          >
            <Group left={margin.left} top={margin.top}>
              {children({ innerWidth, innerHeight, margin })}
            </Group>
          </svg>
        )}
      </div>
      {caption && (
        <p className="mt-1 text-xs text-muted-foreground">{caption}</p>
      )}
    </figure>
  );
}
