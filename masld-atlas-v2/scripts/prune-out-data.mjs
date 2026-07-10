#!/usr/bin/env node
/**
 * Remove `out/data/` from the static export.
 *
 * Next.js (`output: "export"`) copies the ENTIRE `public/` tree — including the
 * ~800 MB `public/data/` (parquets + legacy `genes/`/`network/gene_graphs/`
 * JSONs) — into `out/` at build time. Setting `NEXT_PUBLIC_DATA_BASE` only
 * redirects the browser's *runtime* fetches to the HF Dataset CDN; it does NOT
 * stop the build from bundling the local copy.
 *
 * For the HF **Static Space** deploy, the data lives on the Dataset CDN, so
 * `out/data/` is dead weight that would (a) blow the free-tier size budget and
 * (b) defeat the whole decoupled-data architecture. This script deletes it after
 * a CDN-targeted build. Run ONLY when the build baked a CDN `NEXT_PUBLIC_DATA_BASE`
 * (i.e. via `npm run build:cdn`) — a self-contained `npm run build` needs `out/data`.
 */
import { rm, stat } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const projectRoot = join(dirname(fileURLToPath(import.meta.url)), "..");
const outData = join(projectRoot, "out", "data");

try {
  await stat(outData);
} catch {
  console.log("[prune-out-data] out/data/ not present — nothing to prune.");
  process.exit(0);
}

await rm(outData, { recursive: true, force: true });
console.log(
  "[prune-out-data] Removed out/data/ — the Static Space serves data from the " +
    "HF Dataset CDN (NEXT_PUBLIC_DATA_BASE). Upload out/ to the Space now."
);
