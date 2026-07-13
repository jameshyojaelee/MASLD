"use client";

/**
 * Diverging bar chart: signed values as horizontal bars centered on zero
 * (e.g. per-gene logFC, or a female−male effect contrast). Direction is encoded
 * by the diverging hue (blue = negative, red = positive); significance by
 * opacity. An optional per-row control anchor is drawn as a gray tick.
 */

import { useMemo } from "react";
import { scaleLinear, scaleBand } from "@visx/scale";
import { Line } from "@visx/shape";
import { Text } from "@visx/text";
import { divergingColor, CONTROL, NONSIG_OPACITY } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { BottomAxis, GridCols, AXIS_TEXT, MarkGroup } from "./chart-primitives";

export interface DivergingBarRow {
  label: string;
  /** Signed magnitude; bar extends left (neg) / right (pos) from zero. */
  value: number;
  /** Significant? Controls bar opacity. Defaults to true. */
  sig?: boolean;
  /** Optional control-anchor value, drawn as a gray reference tick. */
  baseline?: number;
  /** Optional explicit color override (else diverging-by-sign). */
  color?: string;
}

export interface DivergingBarProps {
  data: DivergingBarRow[];
  xLabel?: string;
  /** Magnitude that saturates full blue/red in the diverging ramp. */
  clamp?: number;
  valueFormat?: (v: number) => string;
  /** Show the numeric value at each bar's end. */
  showValues?: boolean;
  title?: React.ReactNode;
  caption?: React.ReactNode;
  height?: number;
  ariaLabel?: string;
}

export function DivergingBar({
  data,
  xLabel = "log₂ fold-change",
  clamp = 2,
  valueFormat = (v) => v.toFixed(2),
  showValues = false,
  title,
  caption,
  height,
  ariaLabel,
}: DivergingBarProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<DivergingBarRow>();

  const xDomain = useMemo(() => {
    const vals: number[] = [0];
    for (const d of data) {
      vals.push(d.value);
      if (d.baseline != null) vals.push(d.baseline);
    }
    const maxAbs = Math.max(1e-6, ...vals.map((v) => Math.abs(v)));
    const pad = maxAbs * 0.08;
    return [-(maxAbs + pad), maxAbs + pad] as [number, number];
  }, [data]);

  const resolvedHeight = height ?? Math.max(120, data.length * 26 + 52);
  const rowLabelW = 112;

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={resolvedHeight}
        ariaLabel={ariaLabel ?? "Diverging bar chart"}
        margin={{ top: 8, right: 24, bottom: 40, left: rowLabelW }}
      >
        {({ innerWidth, innerHeight }) => {
          const xScale = scaleLinear({ domain: xDomain, range: [0, innerWidth] });
          const yScale = scaleBand({
            domain: data.map((d) => d.label),
            range: [0, innerHeight],
            padding: 0.28,
          });
          const barH = yScale.bandwidth();
          const x0 = xScale(0);
          const xTicks = xScale.ticks(6);

          return (
            <>
              <GridCols scale={xScale} ticks={xTicks} height={innerHeight} />

              {/* Zero anchor line in control gray */}
              <Line
                from={{ x: x0, y: 0 }}
                to={{ x: x0, y: innerHeight }}
                stroke={CONTROL}
                strokeWidth={1.25}
              />

              <MarkGroup>
              {data.map((d) => {
                const y = yScale(d.label) ?? 0;
                const sig = d.sig ?? true;
                const xv = xScale(d.value);
                const left = Math.min(x0, xv);
                const w = Math.abs(xv - x0);
                const fill = d.color ?? divergingColor(d.value, clamp);
                return (
                  <g
                    key={d.label}
                    onMouseMove={(e) => show(e, d)}
                    onMouseLeave={hide}
                  >
                    <rect
                      x={left}
                      y={y}
                      width={w}
                      height={barH}
                      rx={2}
                      fill={fill}
                      fillOpacity={sig ? 1 : NONSIG_OPACITY}
                    />
                    {/* Control-anchor tick */}
                    {d.baseline != null && (
                      <Line
                        from={{ x: xScale(d.baseline), y: y - 1 }}
                        to={{ x: xScale(d.baseline), y: y + barH + 1 }}
                        stroke={CONTROL}
                        strokeWidth={2}
                      />
                    )}
                    <Text
                      x={-8}
                      y={y + barH / 2}
                      verticalAnchor="middle"
                      textAnchor="end"
                      fontSize={10}
                      fill={AXIS_TEXT}
                      width={rowLabelW - 12}
                    >
                      {d.label}
                    </Text>
                    {showValues && (
                      <Text
                        x={xv + (d.value >= 0 ? 4 : -4)}
                        y={y + barH / 2}
                        verticalAnchor="middle"
                        textAnchor={d.value >= 0 ? "start" : "end"}
                        fontSize={9}
                        fill={AXIS_TEXT}
                      >
                        {valueFormat(d.value)}
                      </Text>
                    )}
                  </g>
                );
              })}
              </MarkGroup>

              <BottomAxis
                scale={xScale}
                top={innerHeight}
                numTicks={6}
                label={xLabel}
              />
            </>
          );
        }}
      </ChartFrame>

      {tooltip && (
        <ChartTooltip left={tooltip.left} top={tooltip.top}>
          <div className="font-medium">{tooltip.data.label}</div>
          <div>{valueFormat(tooltip.data.value)}</div>
          {tooltip.data.baseline != null && (
            <div className="text-muted-foreground">
              control {valueFormat(tooltip.data.baseline)}
            </div>
          )}
        </ChartTooltip>
      )}
    </div>
  );
}
