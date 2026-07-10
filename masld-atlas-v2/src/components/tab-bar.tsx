"use client";

import * as React from "react";
import type { ReactNode } from "react";
import { m, useReducedMotion } from "framer-motion";
import { cn } from "@/lib/utils";

export interface TabItem {
  id: string;
  label: ReactNode;
}

export interface TabBarProps {
  tabs: TabItem[];
  active: string;
  onChange: (id: string) => void;
  className?: string;
  "aria-label"?: string;
}

/**
 * Accessible underline tab bar. Roving-tabindex + arrow/Home/End keyboard
 * navigation (automatic activation), following the WAI-ARIA tabs pattern.
 *
 * The active-tab underline is a single element that slides between tabs. It is
 * measured off the active button (offsetLeft/offsetWidth) and animated via
 * transform + width — this works under the app's `LazyMotion(domAnimation)`
 * feature set (which does NOT ship the `layout`/`layoutId` feature) and snaps
 * instantly under `prefers-reduced-motion` (global MotionConfig + local guard).
 */
export function TabBar({
  tabs,
  active,
  onChange,
  className,
  "aria-label": ariaLabel,
}: TabBarProps) {
  const refs = React.useRef<(HTMLButtonElement | null)[]>([]);
  const listRef = React.useRef<HTMLDivElement>(null);
  const reduce = useReducedMotion();
  const [indicator, setIndicator] = React.useState<{
    left: number;
    width: number;
    ready: boolean;
  }>({ left: 0, width: 0, ready: false });

  const measure = React.useCallback(() => {
    const idx = tabs.findIndex((t) => t.id === active);
    const el = refs.current[idx];
    if (!el || !listRef.current) return;
    setIndicator({ left: el.offsetLeft, width: el.offsetWidth, ready: true });
  }, [tabs, active]);

  React.useLayoutEffect(() => {
    measure();
  }, [measure]);

  React.useEffect(() => {
    const el = listRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver(() => measure());
    ro.observe(el);
    return () => ro.disconnect();
  }, [measure]);

  const focusTab = (index: number) => {
    const next = (index + tabs.length) % tabs.length;
    const tab = tabs[next];
    if (!tab) return;
    onChange(tab.id);
    refs.current[next]?.focus();
  };

  const handleKeyDown = (e: React.KeyboardEvent, index: number) => {
    switch (e.key) {
      case "ArrowRight":
      case "ArrowDown":
        e.preventDefault();
        focusTab(index + 1);
        break;
      case "ArrowLeft":
      case "ArrowUp":
        e.preventDefault();
        focusTab(index - 1);
        break;
      case "Home":
        e.preventDefault();
        focusTab(0);
        break;
      case "End":
        e.preventDefault();
        focusTab(tabs.length - 1);
        break;
    }
  };

  return (
    <div
      ref={listRef}
      role="tablist"
      aria-label={ariaLabel}
      aria-orientation="horizontal"
      className={cn(
        "relative flex items-center gap-1 overflow-x-auto border-b border-border",
        className
      )}
    >
      {tabs.map((tab, i) => {
        const isActive = tab.id === active;
        return (
          <button
            key={tab.id}
            ref={(el) => {
              refs.current[i] = el;
            }}
            type="button"
            role="tab"
            id={`tab-${tab.id}`}
            aria-selected={isActive}
            tabIndex={isActive ? 0 : -1}
            onClick={() => onChange(tab.id)}
            onKeyDown={(e) => handleKeyDown(e, i)}
            className={cn(
              "relative -mb-px shrink-0 rounded-t-sm border-b-2 border-transparent px-3 py-2 text-sm font-medium whitespace-nowrap transition-colors outline-none focus-visible:ring-2 focus-visible:ring-ring/50",
              isActive
                ? "text-foreground"
                : "text-muted-foreground hover:border-border hover:text-foreground"
            )}
          >
            {tab.label}
          </button>
        );
      })}
      <m.span
        aria-hidden
        className="pointer-events-none absolute bottom-0 left-0 h-0.5 rounded-full bg-primary"
        initial={false}
        animate={{ x: indicator.left, width: indicator.width }}
        transition={
          reduce
            ? { duration: 0 }
            : { type: "spring", stiffness: 420, damping: 34, mass: 0.7 }
        }
        style={{ opacity: indicator.ready ? 1 : 0 }}
      />
    </div>
  );
}
