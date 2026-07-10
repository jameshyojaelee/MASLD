"use client";

import dynamic from "next/dynamic";
import type { ComponentType } from "react";
import { useHashLocation } from "@/lib/hash-router";

// Route -> page component. Each page is lazy (ssr:false) so it stays a separate
// client chunk loaded on demand — the landing bundle is not inflated. Only `/`
// (index.html) is ever served by the static host; every other route renders
// here, driven by location.hash. Home ("/") renders `children` (the prerendered
// root page) untouched.
function Loading() {
  return (
    <div className="mx-auto w-full max-w-6xl px-6 py-8">
      <div className="h-8 w-48 animate-pulse rounded bg-muted/50" />
      <div className="mt-6 h-64 animate-pulse rounded-lg bg-muted/40" />
    </div>
  );
}

const ROUTES: Record<string, ComponentType> = {
  "/atlas": dynamic(() => import("@/app/atlas/page"), { ssr: false, loading: Loading }),
  "/explore": dynamic(() => import("@/app/explore/page"), { ssr: false, loading: Loading }),
  "/single-cell": dynamic(() => import("@/app/single-cell/page"), { ssr: false, loading: Loading }),
  "/genetics": dynamic(() => import("@/app/genetics/page"), { ssr: false, loading: Loading }),
  "/progression": dynamic(() => import("@/app/progression/page"), { ssr: false, loading: Loading }),
  "/drugs": dynamic(() => import("@/app/drugs/page"), { ssr: false, loading: Loading }),
  "/programs": dynamic(() => import("@/app/programs/page"), { ssr: false, loading: Loading }),
  "/translation": dynamic(() => import("@/app/translation/page"), { ssr: false, loading: Loading }),
  "/proteomics": dynamic(() => import("@/app/proteomics/page"), { ssr: false, loading: Loading }),
  "/network": dynamic(() => import("@/app/network/page"), { ssr: false, loading: Loading }),
  "/downloads": dynamic(() => import("@/app/downloads/page"), { ssr: false, loading: Loading }),
  "/species": dynamic(() => import("@/app/species/page"), { ssr: false, loading: Loading }),
  "/gene": dynamic(() => import("@/app/gene/page"), { ssr: false, loading: Loading }),
  "/causal": dynamic(() => import("@/app/causal/page"), { ssr: false, loading: Loading }),
  "/resolution": dynamic(() => import("@/app/resolution/page"), { ssr: false, loading: Loading }),
};

export function RouteSwitch({ children }: { children: React.ReactNode }) {
  const { path } = useHashLocation();
  if (path === "/") return <>{children}</>;
  const View = ROUTES[path];
  return View ? <View /> : <>{children}</>;
}
