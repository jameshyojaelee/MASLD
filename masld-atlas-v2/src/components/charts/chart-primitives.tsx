"use client";

/**
 * Shared visual primitives for the chart library: themed axes, subtle
 * gridlines, and a small legend row. Centralising these keeps every chart on
 * one visual system (axis weight, tick size, muted label color) so the portal
 * reads as a single scientific instrument.
 *
 * All text uses theme foreground/muted tokens (never a chart hue) per project
 * doctrine: figure text is never colored.
 */

import { AxisBottom, AxisLeft } from "@visx/axis";
import type { AxisScale } from "@visx/axis";
import { type ReactNode } from "react";
import { m, useReducedMotion } from "framer-motion";
import type { Variants } from "framer-motion";

// Theme tokens resolved as CSS-var strings. These render in SVG (not canvas),
// so CSS custom properties are honored and adapt to light/dark automatically.
export const AXIS_LINE = "var(--color-border)";
export const AXIS_TEXT = "var(--color-muted-foreground)";
export const GRID_LINE = "var(--color-border)";

export const TICK_FONT_SIZE = 10;
export const LABEL_FONT_SIZE = 11;

/** Standard tick-label props for the horizontal (bottom) axis. */
export const bottomTickLabelProps = {
  fill: AXIS_TEXT,
  fontSize: TICK_FONT_SIZE,
  fontFamily: "inherit",
  textAnchor: "middle" as const,
  dy: "0.25em",
};

/** Standard tick-label props for the vertical (left) axis. */
export const leftTickLabelProps = {
  fill: AXIS_TEXT,
  fontSize: TICK_FONT_SIZE,
  fontFamily: "inherit",
  textAnchor: "end" as const,
  dx: "-0.25em",
  dy: "0.25em",
};

const axisLabelProps = {
  fill: AXIS_TEXT,
  fontSize: LABEL_FONT_SIZE,
  fontFamily: "inherit",
  textAnchor: "middle" as const,
};

// ---------------------------------------------------------------------------
// Gridlines
// ---------------------------------------------------------------------------

/** Horizontal gridlines at the supplied y tick values. */
export function GridRows({ scale, ticks, width }: {
  scale: (v: number) => number | undefined;
  ticks: number[];
  width: number;
}) {
  return (
    <g aria-hidden>
      {ticks.map((t, i) => {
        const y = scale(t);
        if (y == null) return null;
        return (
          <line
            key={i}
            x1={0}
            x2={width}
            y1={y}
            y2={y}
            stroke={GRID_LINE}
            strokeOpacity={0.35}
            strokeWidth={1}
            shapeRendering="crispEdges"
          />
        );
      })}
    </g>
  );
}

/** Vertical gridlines at the supplied x tick values. */
export function GridCols({ scale, ticks, height }: {
  scale: (v: number) => number | undefined;
  ticks: number[];
  height: number;
}) {
  return (
    <g aria-hidden>
      {ticks.map((t, i) => {
        const x = scale(t);
        if (x == null) return null;
        return (
          <line
            key={i}
            x1={x}
            x2={x}
            y1={0}
            y2={height}
            stroke={GRID_LINE}
            strokeOpacity={0.35}
            strokeWidth={1}
            shapeRendering="crispEdges"
          />
        );
      })}
    </g>
  );
}

// ---------------------------------------------------------------------------
// Themed axes (thin wrappers over visx with the shared styling baked in)
// ---------------------------------------------------------------------------

export function BottomAxis({
  scale,
  top,
  label,
  numTicks,
  tickFormat,
  hideAxisLine,
}: {
  scale: AxisScale;
  top: number;
  label?: string;
  numTicks?: number;
  tickFormat?: (v: never, i: number) => string;
  hideAxisLine?: boolean;
}) {
  return (
    <AxisBottom
      scale={scale}
      top={top}
      stroke={AXIS_LINE}
      tickStroke={AXIS_LINE}
      hideAxisLine={hideAxisLine}
      numTicks={numTicks}
      tickFormat={tickFormat as never}
      tickLabelProps={() => bottomTickLabelProps}
      label={label}
      labelProps={axisLabelProps}
      labelOffset={20}
    />
  );
}

export function LeftAxis({
  scale,
  label,
  numTicks,
  tickFormat,
  hideAxisLine,
}: {
  scale: AxisScale;
  label?: string;
  numTicks?: number;
  tickFormat?: (v: never, i: number) => string;
  hideAxisLine?: boolean;
}) {
  return (
    <AxisLeft
      scale={scale}
      stroke={AXIS_LINE}
      tickStroke={AXIS_LINE}
      hideAxisLine={hideAxisLine}
      numTicks={numTicks}
      tickFormat={tickFormat as never}
      tickLabelProps={() => leftTickLabelProps}
      label={label}
      labelProps={{ ...axisLabelProps, dy: "-1.2em" }}
      labelOffset={28}
    />
  );
}

// ---------------------------------------------------------------------------
// Legend
// ---------------------------------------------------------------------------

export interface LegendItem {
  label: string;
  color: string;
  /** Optional shape override; default is a filled swatch. */
  shape?: "swatch" | "line";
}

/**
 * Compact wrapping legend row rendered in the DOM (below the SVG).
 *
 * Backward-compatible: with only `items`/`className` it renders static swatch
 * spans exactly as before. Passing `onItemClick` promotes each entry to a
 * click-to-isolate toggle button; entries in `mutedLabels` render dimmed (their
 * series hidden by the consumer). Interaction is opt-in and non-breaking.
 */
export function ChartLegend({
  items,
  className,
  mutedLabels,
  onItemClick,
}: {
  items: LegendItem[];
  className?: string;
  /** Labels rendered dimmed (their series toggled off by the consumer). */
  mutedLabels?: Set<string>;
  /** When set, entries become click-to-isolate toggle buttons. */
  onItemClick?: (label: string, index: number) => void;
}) {
  const swatch = (it: LegendItem) =>
    it.shape === "line" ? (
      <span
        className="inline-block h-0.5 w-3.5 rounded-full"
        style={{ backgroundColor: it.color }}
      />
    ) : (
      <span
        className="inline-block h-2.5 w-2.5 rounded-[2px]"
        style={{ backgroundColor: it.color }}
      />
    );

  return (
    <div
      className={
        "mt-2 flex flex-wrap items-center gap-x-3 gap-y-1 text-xs text-muted-foreground " +
        (className ?? "")
      }
    >
      {items.map((it, i) => {
        const muted = mutedLabels?.has(it.label) ?? false;
        if (onItemClick) {
          return (
            <button
              key={i}
              type="button"
              onClick={() => onItemClick(it.label, i)}
              aria-pressed={!muted}
              className={
                "inline-flex items-center gap-1.5 rounded-sm px-1 py-0.5 transition-opacity hover:bg-muted/60 " +
                (muted ? "opacity-40" : "opacity-100")
              }
            >
              {swatch(it)}
              <span className="tabular-nums">{it.label}</span>
            </button>
          );
        }
        return (
          <span key={i} className="inline-flex items-center gap-1.5">
            {swatch(it)}
            <span className="tabular-nums">{it.label}</span>
          </span>
        );
      })}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Misc helpers
// ---------------------------------------------------------------------------

/** Wrap a chart body so an absolutely-positioned tooltip anchors correctly. */
export function ChartCanvas({
  wrapperRef,
  children,
}: {
  wrapperRef: React.RefObject<HTMLDivElement | null>;
  children: ReactNode;
}) {
  return (
    <div ref={wrapperRef} className="relative w-full">
      {children}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Mark-mount animation
// ---------------------------------------------------------------------------

/**
 * Above this many marks, per-mark staggered entry is disabled and the whole
 * group cross-fades instead — keeps large scatters (thousands of genes) from
 * spawning thousands of animation nodes on the free-tier bundle.
 */
export const MARK_STAGGER_MAX = 400;

/**
 * Per-mark entry variant (opacity only). CRITICAL: this animates the *element*
 * opacity 0→1 (the approach), never the mark's resting `fillOpacity`, which is
 * the significance encoding and must stay put. `custom` = mark index for the
 * stagger delay. Consume as `<m.circle custom={i} variants={markItemVariants}
 * initial="hidden" animate="show" fillOpacity={sigOpacity} />`.
 */
export const markItemVariants: Variants = {
  hidden: { opacity: 0 },
  show: (i: number) => ({
    opacity: 1,
    transition: {
      delay: Math.min(i * 0.003, 0.6),
      duration: 0.25,
      ease: [0.16, 1, 0.3, 1],
    },
  }),
};

/**
 * Wraps a chart's marks in a single group cross-fade on first mount. O(1) in
 * mark count (unlike per-mark stagger), so it is the safe default for large or
 * unbounded mark sets. Renders a plain <g> under reduced motion, so the resting
 * frame is always the final (correctly-encoded) frame.
 */
export function MarkGroup({
  children,
  className,
}: {
  children: ReactNode;
  className?: string;
}) {
  const reduce = useReducedMotion();
  if (reduce) return <g className={className}>{children}</g>;
  return (
    <m.g
      className={className}
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.35, ease: [0.16, 1, 0.3, 1] }}
    >
      {children}
    </m.g>
  );
}
