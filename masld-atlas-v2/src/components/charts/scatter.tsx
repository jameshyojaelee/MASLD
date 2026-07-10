"use client";

/**
 * Categorical scatter: two continuous axes with optional x=0 / y=0 reference
 * lines and quadrant annotations.
 *
 * Every mark color must be supplied already palette-derived (via
 * `@/lib/palette` — e.g. `speciesCategoryColor`, `categoricalColor`), with
 * control / non-significant categories resolving to `CONTROL`. Inactive marks
 * dim by opacity (never by a different hue). Axis + quadrant text uses theme
 * tokens only, never a data hue.
 */

import { useMemo, type ReactNode } from "react";
import { scaleLinear } from "@visx/scale";
import { Line } from "@visx/shape";
import { Text } from "@visx/text";
import { CONTROL } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import {
  BottomAxis,
  LeftAxis,
  GridRows,
  GridCols,
  AXIS_LINE,
  AXIS_TEXT,
} from "./chart-primitives";

export interface ScatterPoint {
  x: number;
  y: number;
  label: string;
  /** Palette-derived mark color. Defaults to `CONTROL`. */
  color?: string;
  /** Free-form category (for the caller's tooltip / grouping). */
  category?: string;
  /** Render dimmed (inactive filter). */
  dimmed?: boolean;
  /** Render with the emphasized (larger) radius. */
  emphasized?: boolean;
}

export interface QuadrantLabel {
  x: "left" | "right";
  y: "top" | "bottom";
  text: string;
  /** Render at lower emphasis (e.g. the discordant quadrants). */
  muted?: boolean;
}

export interface ScatterProps {
  data: ScatterPoint[];
  xLabel?: string;
  yLabel?: string;
  /** Draw x=0 / y=0 reference lines (default true). */
  refLines?: boolean;
  quadrantLabels?: QuadrantLabel[];
  /** Fractional padding added to the symmetric per-axis domain. */
  domainPad?: number;
  height?: number;
  title?: ReactNode;
  caption?: ReactNode;
  ariaLabel?: string;
  onPointClick?: (label: string) => void;
  /** Tooltip body builder; defaults to label + coordinates. */
  tooltipLines?: (p: ScatterPoint) => ReactNode;
  baseRadius?: number;
  emphasizedRadius?: number;
}

function symmetricDomain(vals: number[], pad: number): [number, number] {
  const maxAbs = Math.max(0.1, ...vals.map((v) => Math.abs(v)));
  const ext = maxAbs * (1 + pad);
  return [-ext, ext];
}

export function Scatter({
  data,
  xLabel,
  yLabel,
  refLines = true,
  quadrantLabels,
  domainPad = 0.08,
  height = 420,
  title,
  caption,
  ariaLabel,
  onPointClick,
  tooltipLines,
  baseRadius = 2.5,
  emphasizedRadius = 3.5,
}: ScatterProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<ScatterPoint>();

  const xDomain = useMemo(
    () => symmetricDomain(data.map((d) => d.x), domainPad),
    [data, domainPad]
  );
  const yDomain = useMemo(
    () => symmetricDomain(data.map((d) => d.y), domainPad),
    [data, domainPad]
  );

  // Draw dimmed marks first so active marks sit on top.
  const ordered = useMemo(
    () => [...data].sort((a, b) => Number(!!b.dimmed) - Number(!!a.dimmed)),
    [data]
  );

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Scatter plot"}
        margin={{ top: 14, right: 16, bottom: 46, left: 54 }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleLinear({ domain: xDomain, range: [0, innerWidth], nice: true });
          const yScale = scaleLinear({ domain: yDomain, range: [innerHeight, 0], nice: true });
          const xTicks = xScale.ticks(6);
          const yTicks = yScale.ticks(6);
          const x0 = xScale(0);
          const y0 = yScale(0);

          return (
            <>
              <GridRows scale={yScale} ticks={yTicks} width={innerWidth} />
              <GridCols scale={xScale} ticks={xTicks} height={innerHeight} />

              {refLines && (
                <>
                  <Line
                    from={{ x: x0, y: 0 }}
                    to={{ x: x0, y: innerHeight }}
                    stroke={AXIS_LINE}
                    strokeWidth={1.25}
                  />
                  <Line
                    from={{ x: 0, y: y0 }}
                    to={{ x: innerWidth, y: y0 }}
                    stroke={AXIS_LINE}
                    strokeWidth={1.25}
                  />
                </>
              )}

              {quadrantLabels?.map((q, i) => (
                <Text
                  key={`q-${i}`}
                  x={innerWidth * (q.x === "left" ? 0.22 : 0.78)}
                  y={innerHeight * (q.y === "top" ? 0.07 : 0.95)}
                  fontSize={10}
                  textAnchor="middle"
                  fill={AXIS_TEXT}
                  opacity={q.muted ? 0.55 : 1}
                >
                  {q.text}
                </Text>
              ))}

              {ordered.map((p, i) => {
                const cx = xScale(p.x);
                const cy = yScale(p.y);
                const r = p.emphasized ? emphasizedRadius : baseRadius;
                return (
                  <circle
                    key={`${p.label}-${i}`}
                    cx={cx}
                    cy={cy}
                    r={p.dimmed ? Math.max(1.5, r - 1) : r}
                    fill={p.color ?? CONTROL}
                    fillOpacity={p.dimmed ? 0.12 : 0.8}
                    className={onPointClick && !p.dimmed ? "cursor-pointer" : undefined}
                    onMouseMove={p.dimmed ? undefined : (e) => show(e, p)}
                    onMouseLeave={p.dimmed ? undefined : hide}
                    onClick={
                      onPointClick && !p.dimmed
                        ? () => onPointClick(p.label)
                        : undefined
                    }
                  />
                );
              })}

              <BottomAxis
                scale={xScale}
                top={innerHeight}
                numTicks={6}
                label={xLabel}
              />
              <LeftAxis scale={yScale} numTicks={6} label={yLabel} />
            </>
          );
        }}
      </ChartFrame>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          {tooltipLines ? (
            tooltipLines(tooltip.data)
          ) : (
            <>
              <div className="font-medium italic">{tooltip.data.label}</div>
              <div>
                x {tooltip.data.x.toFixed(2)} · y {tooltip.data.y.toFixed(2)}
              </div>
            </>
          )}
        </ChartTooltip>
      )}
    </div>
  );
}
