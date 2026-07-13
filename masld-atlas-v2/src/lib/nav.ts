import {
  Home,
  LayoutGrid,
  Compass,
  Microscope,
  Dna,
  TrendingUp,
  Pill,
  Boxes,
  GitMerge,
  FlaskConical,
  Network,
  Download,
  PawPrint,
} from "lucide-react";

/** A lucide icon component (accepts className / size). */
export type NavIcon = React.ComponentType<{
  className?: string;
  size?: number | string;
  "aria-hidden"?: boolean;
}>;

export interface NavItem {
  label: string;
  href: string;
  icon: NavIcon;
  /** Flagged for keep/retire review — rendered with a small marker. */
  flag?: boolean;
}

export interface NavSection {
  title: string;
  items: NavItem[];
}

/**
 * Single source of truth for the primary navigation / information architecture.
 * Both the desktop sidebar and the mobile drawer render from this array, so IA
 * changes happen in one place.
 *
 * NOTE: `/single-cell` will absorb `/resolution` and `/genetics` will absorb
 * `/causal`; those legacy routes remain live for now (migrated later).
 */
export const NAV: NavSection[] = [
  {
    title: "OVERVIEW",
    items: [
      { label: "Home", href: "/", icon: Home },
      { label: "Atlas", href: "/atlas", icon: LayoutGrid },
      { label: "Explore", href: "/explore", icon: Compass },
    ],
  },
  {
    title: "MULTI-SCALE",
    items: [{ label: "Single-Cell", href: "/single-cell", icon: Microscope }],
  },
  {
    title: "GENETICS & ANCESTRY",
    items: [{ label: "Genetics", href: "/genetics", icon: Dna }],
  },
  {
    title: "PROGRESSION",
    items: [{ label: "Progression", href: "/progression", icon: TrendingUp }],
  },
  {
    title: "THERAPEUTICS",
    items: [
      { label: "Drugs", href: "/drugs", icon: Pill },
      { label: "Programs", href: "/programs", icon: Boxes },
    ],
  },
  {
    title: "INTEGRATE",
    items: [
      { label: "Translation", href: "/translation", icon: GitMerge },
      { label: "Proteomics", href: "/proteomics", icon: FlaskConical },
      { label: "Network", href: "/network", icon: Network },
    ],
  },
  {
    title: "TOOLS",
    items: [
      { label: "Downloads", href: "/downloads", icon: Download },
      { label: "Cross-Species", href: "/species", icon: PawPrint, flag: true },
    ],
  },
];

/** True when `href` should be highlighted for the current pathname. */
export function isNavActive(pathname: string, href: string): boolean {
  if (href === "/") return pathname === "/";
  return pathname === href || pathname.startsWith(`${href}/`);
}
