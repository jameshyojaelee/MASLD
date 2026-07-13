"use client";

/**
 * Unsigned magnitude bars (horizontal or vertical).
 *
 * For counts / magnitudes that carry no direction. Significance is encoded by
 * opacity (significant = full, non-significant = `NONSIG_OPACITY`), never by a
 * separate hue. Mark color must be supplied already palette-derived (via
 * `@/lib/palette`); the default is `CONTROL`. Axis + value text uses theme
 * tokens only — never a data hue. (Signed effects belong in `DivergingBar`.)
 */

import { useMemo, type ReactNode } from "react";
import { scaleBand, scaleLinear } from "@visx/scale";
import { Text } from "@visx/text";
import { CONTROL, SIG_OPACITY, NONSIG_OPACITY } from "@/lib/palette";
import { ChartFrame } from "./chart-frame";
import { useChartTooltip, ChartTooltip } from "./chart-tooltip";
import { BottomAxis, LeftAxis, GridRows, GridCols, AXIS_TEXT } from "./chart-primitives";

export interface BarDatum {
  label: string;
  value: number;
  /** Palette-derived color. Defaults to `CONTROL`. */
  color?: string;
  /** Significant → full opacity; else `NONSIG_OPACITY`. Default true. */
  sig?: boolean;
  /** Optional short annotation drawn just past the bar end. */
  annotation?: string;
}

export interface BarProps {
  data: BarDatum[];
  orientation?: "horizontal" | "vertical";
  /** Axis label for the value dimension. */
  valueLabel?: string;
  height?: number;
  title?: ReactNode;
  caption?: ReactNode;
  ariaLabel?: string;
  onBarClick?: (label: string) => void;
  tooltipLines?: (d: BarDatum) => ReactNode;
  valueFormat?: (v: number) => string;
  /** Fixed value-axis maximum; defaults to data max. */
  maxValue?: number;
}

export function Bar({
  data,
  orientation = "horizontal",
  valueLabel,
  height = 300,
  title,
  caption,
  ariaLabel,
  onBarClick,
  tooltipLines,
  valueFormat = (v) => v.toLocaleString(),
  maxValue,
}: BarProps) {
  const { wrapperRef, tooltip, show, hide } = useChartTooltip<BarDatum>();

  const vMax = useMemo(
    () => maxValue ?? Math.max(1, ...data.map((d) => d.value)),
    [data, maxValue]
  );
  const labels = useMemo(() => data.map((d) => d.label), [data]);

  const horizontal = orientation === "horizontal";

  return (
    <div ref={wrapperRef} className="relative w-full">
      <ChartFrame
        title={title}
        caption={caption}
        height={height}
        ariaLabel={ariaLabel ?? "Bar chart"}
        margin={
          horizontal
            ? { top: 8, right: 44, bottom: 40, left: 128 }
            : { top: 8, right: 12, bottom: 52, left: 48 }
        }
      >
        {({ innerWidth, innerHeight }) => {
          const band = scaleBand({
            domain: labels,
            range: horizontal ? [0, innerHeight] : [0, innerWidth],
            padding: 0.28,
          });
          const value = scaleLinear({
            domain: [0, vMax],
            range: horizontal ? [0, innerWidth] : [innerHeight, 0],
            nice: true,
          });
          const bw = band.bandwidth();

          return (
            <>
              {horizontal ? (
                <GridCols scale={value} ticks={value.ticks(5)} height={innerHeight} />
              ) : (
                <GridRows scale={value} ticks={value.ticks(5)} width={innerWidth} />
              )}

              {data.map((d) => {
                const sig = d.sig ?? true;
                const bandPos = band(d.label) ?? 0;
                const fill = d.color ?? CONTROL;
                const opacity = sig ? SIG_OPACITY : NONSIG_OPACITY;

                if (horizontal) {
                  const w = value(d.value);
                  const cy = bandPos + bw / 2;
                  return (
                    <g key={d.label}>
                      <rect
                        x={0}
                        y={bandPos}
                        width={Math.max(0, w)}
                        height={bw}
                        rx={2}
                        fill={fill}
                        fillOpacity={opacity}
                        className={onBarClick ? "cursor-pointer" : undefined}
                        onMouseMove={(e) => show(e, d)}
                        onMouseLeave={hide}
                        onClick={onBarClick ? () => onBarClick(d.label) : undefined}
                      />
                      <Text
                        x={w + 6}
                        y={cy}
                        fontSize={10}
                        verticalAnchor="middle"
                        fill={AXIS_TEXT}
                      >
                        {d.annotation ?? valueFormat(d.value)}
                      </Text>
                    </g>
                  );
                }

                const yTop = value(d.value);
                const cx = bandPos + bw / 2;
                return (
                  <g key={d.label}>
                    <rect
                      x={bandPos}
                      y={yTop}
                      width={bw}
                      height={Math.max(0, innerHeight - yTop)}
                      rx={2}
                      fill={fill}
                      fillOpacity={opacity}
                      className={onBarClick ? "cursor-pointer" : undefined}
                      onMouseMove={(e) => show(e, d)}
                      onMouseLeave={hide}
                      onClick={onBarClick ? () => onBarClick(d.label) : undefined}
                    />
                    <Text
                      x={cx}
                      y={yTop - 4}
                      fontSize={10}
                      textAnchor="middle"
                      fill={AXIS_TEXT}
                    >
                      {d.annotation ?? valueFormat(d.value)}
                    </Text>
                  </g>
                );
              })}

              {horizontal ? (
                <>
                  <LeftAxis scale={band} hideAxisLine />
                  <BottomAxis scale={value} top={innerHeight} numTicks={5} label={valueLabel} />
                </>
              ) : (
                <>
                  <LeftAxis scale={value} numTicks={5} label={valueLabel} />
                  <BottomAxis scale={band} top={innerHeight} hideAxisLine />
                </>
              )}
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
              <div className="font-medium">{tooltip.data.label}</div>
              <div>{valueFormat(tooltip.data.value)}</div>
            </>
          )}
        </ChartTooltip>
      )}
    </div>
  );
}
