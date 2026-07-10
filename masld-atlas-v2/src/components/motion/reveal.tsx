"use client";

/**
 * Scroll-reveal wrappers.
 *
 * <Reveal> fades + rises its children into view once (staggering any
 * <RevealItem> children). Honors reduced motion globally via MotionConfig, so
 * no per-instance guard is needed. Built on the lightweight `m` component so it
 * stays inside the LazyMotion feature bundle.
 */

import { m } from "framer-motion";
import type { ReactNode } from "react";
import { staggerContainer, staggerItem, fadeInUp } from "@/lib/motion";

export function Reveal({
  children,
  className,
  stagger = true,
  once = true,
}: {
  children: ReactNode;
  className?: string;
  /** Stagger direct <RevealItem> children (default) vs. a single fade. */
  stagger?: boolean;
  once?: boolean;
}) {
  return (
    <m.div
      className={className}
      variants={stagger ? staggerContainer : fadeInUp}
      initial="hidden"
      whileInView="show"
      viewport={{ once, margin: "-60px" }}
    >
      {children}
    </m.div>
  );
}

export function RevealItem({
  children,
  className,
}: {
  children: ReactNode;
  className?: string;
}) {
  return (
    <m.div className={className} variants={staggerItem}>
      {children}
    </m.div>
  );
}
