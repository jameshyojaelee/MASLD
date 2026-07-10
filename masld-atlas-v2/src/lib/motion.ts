/**
 * Shared motion vocabulary for the portal.
 *
 * Pure data (variants + transition presets) — safe to import from server or
 * client modules. Motion is opt-in and MUST degrade under reduced motion; the
 * global <MotionProvider> sets `MotionConfig reducedMotion="user"`, and
 * per-effect guards live in the components that consume these.
 *
 * Timing mirrors the CSS motion tokens in globals.css so hand-written CSS
 * (`.hover-lift`) and framer animations feel identical.
 */

import type { Variants, Transition } from "framer-motion";

/** Easing + duration constants (seconds), mirroring the CSS tokens. */
export const MOTION = {
  easeOutExpo: [0.16, 1, 0.3, 1] as [number, number, number, number],
  easeSpring: [0.34, 1.56, 0.64, 1] as [number, number, number, number],
  fast: 0.15,
  base: 0.3,
  slow: 0.6,
} as const;

export const springSoft: Transition = {
  type: "spring",
  stiffness: 220,
  damping: 26,
  mass: 0.9,
};

export const easeOut: Transition = {
  duration: MOTION.base,
  ease: MOTION.easeOutExpo,
};

/** Fade + rise. Use as a stagger child or standalone. */
export const fadeInUp: Variants = {
  hidden: { opacity: 0, y: 12 },
  show: { opacity: 1, y: 0, transition: easeOut },
};

/** Container that staggers its children's `show` state. */
export const staggerContainer: Variants = {
  hidden: {},
  show: {
    transition: { staggerChildren: 0.06, delayChildren: 0.02 },
  },
};

/** Child of `staggerContainer`. */
export const staggerItem: Variants = fadeInUp;

/** Small cross-fade for tab-panel content. */
export const tabFade: Variants = {
  hidden: { opacity: 0, y: 4 },
  show: { opacity: 1, y: 0, transition: { duration: MOTION.fast, ease: MOTION.easeOutExpo } },
};

/** Skeleton → content cross-fade. */
export const skeletonToContent: Variants = {
  hidden: { opacity: 0 },
  show: { opacity: 1, transition: { duration: MOTION.base, ease: MOTION.easeOutExpo } },
};
