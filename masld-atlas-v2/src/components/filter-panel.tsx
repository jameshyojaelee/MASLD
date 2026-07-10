"use client";

import * as React from "react";
import type { ReactNode } from "react";
import { ChevronDown } from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";

// ---------------------------------------------------------------------------
// Collapsible section
// ---------------------------------------------------------------------------

export function FilterSection({
  label,
  children,
  defaultOpen = true,
  count,
  className,
}: {
  label: ReactNode;
  children: ReactNode;
  defaultOpen?: boolean;
  /** Optional badge count shown next to the label. */
  count?: number;
  className?: string;
}) {
  const [open, setOpen] = React.useState(defaultOpen);
  return (
    <div className={cn("border-b border-border py-2 last:border-0", className)}>
      <button
        type="button"
        aria-expanded={open}
        onClick={() => setOpen((o) => !o)}
        className="flex w-full items-center justify-between gap-2 py-1 text-left text-xs font-semibold uppercase tracking-wide text-muted-foreground transition-colors hover:text-foreground"
      >
        <span>
          {label}
          {count != null && (
            <span className="ml-1 tabular-nums text-muted-foreground/70">
              ({count})
            </span>
          )}
        </span>
        <ChevronDown
          className={cn(
            "size-3.5 shrink-0 transition-transform",
            open ? "" : "-rotate-90"
          )}
        />
      </button>
      {open && <div className="mt-2 space-y-1.5">{children}</div>}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Dual-thumb range slider (native inputs → static-export-safe, no deps)
// ---------------------------------------------------------------------------

const THUMB =
  "pointer-events-none absolute inset-0 h-4 w-full appearance-none bg-transparent " +
  "[&::-webkit-slider-thumb]:pointer-events-auto [&::-webkit-slider-thumb]:size-3.5 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:border [&::-webkit-slider-thumb]:border-primary [&::-webkit-slider-thumb]:bg-background [&::-webkit-slider-thumb]:shadow-sm " +
  "[&::-moz-range-thumb]:pointer-events-auto [&::-moz-range-thumb]:size-3.5 [&::-moz-range-thumb]:appearance-none [&::-moz-range-thumb]:rounded-full [&::-moz-range-thumb]:border [&::-moz-range-thumb]:border-primary [&::-moz-range-thumb]:bg-background";

export function RangeSlider({
  min,
  max,
  step = 1,
  value,
  onChange,
  unit = "",
  className,
}: {
  min: number;
  max: number;
  step?: number;
  value: [number, number];
  onChange: (value: [number, number]) => void;
  unit?: string;
  className?: string;
}) {
  const [lo, hi] = value;
  const span = max - min || 1;
  const pct = (v: number) => ((v - min) / span) * 100;

  return (
    <div className={cn("w-full", className)}>
      <div className="relative h-4">
        <div className="absolute top-1/2 h-1 w-full -translate-y-1/2 rounded-full bg-muted" />
        <div
          className="absolute top-1/2 h-1 -translate-y-1/2 rounded-full bg-primary"
          style={{ left: `${pct(lo)}%`, right: `${100 - pct(hi)}%` }}
        />
        <input
          type="range"
          aria-label="Minimum"
          min={min}
          max={max}
          step={step}
          value={lo}
          onChange={(e) =>
            onChange([Math.min(Number(e.target.value), hi), hi])
          }
          className={THUMB}
        />
        <input
          type="range"
          aria-label="Maximum"
          min={min}
          max={max}
          step={step}
          value={hi}
          onChange={(e) =>
            onChange([lo, Math.max(Number(e.target.value), lo)])
          }
          className={THUMB}
        />
      </div>
      <div className="mt-1 flex justify-between text-[11px] tabular-nums text-muted-foreground">
        <span>
          {lo}
          {unit}
        </span>
        <span>
          {hi}
          {unit}
        </span>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Schema-driven FilterPanel
// ---------------------------------------------------------------------------

export interface FilterOption {
  value: string;
  label: ReactNode;
  count?: number;
}

export type FilterControl =
  | {
      type: "checkbox";
      key: string;
      label: string;
      options: FilterOption[];
      defaultOpen?: boolean;
    }
  | {
      type: "select";
      key: string;
      label: string;
      options: FilterOption[];
      placeholder?: string;
      defaultOpen?: boolean;
    }
  | {
      type: "range";
      key: string;
      label: string;
      min: number;
      max: number;
      step?: number;
      unit?: string;
      defaultOpen?: boolean;
    };

export type FilterValue = string[] | string | [number, number];

function defaultValueFor(control: FilterControl): FilterValue {
  if (control.type === "checkbox") return [];
  if (control.type === "select") return "";
  return [control.min, control.max];
}

export interface FilterPanelProps {
  controls: FilterControl[];
  values: Record<string, FilterValue | undefined>;
  onChange: (key: string, value: FilterValue | undefined) => void;
  /** Called by "Clear all". If omitted, each control is reset to its empty value. */
  onClearAll?: () => void;
  title?: ReactNode;
  className?: string;
}

/**
 * Faceted filter panel: collapsible checkbox groups, select dropdowns, and a
 * dual-thumb range slider, fully controlled via `values` / `onChange`.
 */
export function FilterPanel({
  controls,
  values,
  onChange,
  onClearAll,
  title = "Filters",
  className,
}: FilterPanelProps) {
  const handleClearAll = () => {
    if (onClearAll) {
      onClearAll();
      return;
    }
    for (const control of controls) {
      onChange(control.key, defaultValueFor(control));
    }
  };

  return (
    <div
      className={cn(
        "rounded-lg border border-border bg-card px-4 py-3 shadow-sm",
        className
      )}
    >
      <div className="flex items-center justify-between gap-2 pb-1">
        <span className="text-sm font-semibold">{title}</span>
        <Button variant="ghost" size="xs" onClick={handleClearAll}>
          Clear all
        </Button>
      </div>

      {controls.map((control) => (
        <FilterSection
          key={control.key}
          label={control.label}
          defaultOpen={control.defaultOpen ?? true}
        >
          {control.type === "checkbox" && (
            <CheckboxGroup
              options={control.options}
              value={(values[control.key] as string[]) ?? []}
              onChange={(v) => onChange(control.key, v)}
            />
          )}

          {control.type === "select" && (
            <select
              value={(values[control.key] as string) ?? ""}
              onChange={(e) => onChange(control.key, e.target.value)}
              className="h-8 w-full rounded-lg border border-input bg-transparent px-2.5 text-sm outline-none transition-colors focus-visible:border-ring focus-visible:ring-3 focus-visible:ring-ring/50 dark:bg-input/30"
            >
              <option value="">{control.placeholder ?? "All"}</option>
              {control.options.map((opt) => (
                <option key={opt.value} value={opt.value}>
                  {typeof opt.label === "string" ? opt.label : opt.value}
                </option>
              ))}
            </select>
          )}

          {control.type === "range" && (
            <RangeSlider
              min={control.min}
              max={control.max}
              step={control.step}
              unit={control.unit}
              value={
                (values[control.key] as [number, number]) ?? [
                  control.min,
                  control.max,
                ]
              }
              onChange={(v) => onChange(control.key, v)}
            />
          )}
        </FilterSection>
      ))}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Checkbox group
// ---------------------------------------------------------------------------

export function CheckboxGroup({
  options,
  value,
  onChange,
}: {
  options: FilterOption[];
  value: string[];
  onChange: (value: string[]) => void;
}) {
  const toggle = (v: string) => {
    onChange(value.includes(v) ? value.filter((x) => x !== v) : [...value, v]);
  };
  return (
    <div className="space-y-0.5">
      {options.map((opt) => (
        <label
          key={opt.value}
          className="flex cursor-pointer items-center gap-2 rounded-md px-1 py-1 text-sm transition-colors hover:bg-muted/50"
        >
          <input
            type="checkbox"
            checked={value.includes(opt.value)}
            onChange={() => toggle(opt.value)}
            className="size-3.5 shrink-0 rounded accent-primary"
          />
          <span className="flex-1 truncate">{opt.label}</span>
          {opt.count != null && (
            <span className="tabular-nums text-xs text-muted-foreground">
              {opt.count}
            </span>
          )}
        </label>
      ))}
    </div>
  );
}
