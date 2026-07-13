#!/usr/bin/env node
/**
 * Build guard: fail if the banned method name appears anywhere in `src/`.
 *
 * Project doctrine forbids the word in code, labels, and copy. The canonical
 * DEG gate is described only as the "effect-size-aware interval-null FDR gate
 * (fdr < 0.05 at lfc = 0.25)". This runs as `prebuild` and `pretest` so a
 * regression can never ship.
 *
 * Match rule: case-insensitive, word-boundary — flags the standalone word
 * (TREAT / Treat / treat) but NOT "treatment" / "treated".
 */
import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, extname, relative } from "node:path";

const ROOT = process.cwd();
const SRC = join(ROOT, "src");
const EXTS = new Set([
  ".ts",
  ".tsx",
  ".js",
  ".jsx",
  ".mjs",
  ".cjs",
  ".css",
  ".md",
  ".mdx",
  ".json",
]);
const BANNED = /\btreat\b/i;

/** @param {string} dir @param {string[]} out */
function walk(dir, out) {
  for (const name of readdirSync(dir)) {
    const full = join(dir, name);
    const st = statSync(full);
    if (st.isDirectory()) {
      if (name === "node_modules" || name.startsWith(".")) continue;
      walk(full, out);
    } else if (EXTS.has(extname(name))) {
      out.push(full);
    }
  }
}

const files = [];
try {
  walk(SRC, files);
} catch (e) {
  console.error(`[check-no-banned-term] cannot scan ${SRC}: ${e.message}`);
  process.exit(2);
}

const hits = [];
for (const file of files) {
  const lines = readFileSync(file, "utf8").split(/\r?\n/);
  lines.forEach((line, i) => {
    if (BANNED.test(line)) {
      hits.push(`${relative(ROOT, file)}:${i + 1}: ${line.trim()}`);
    }
  });
}

if (hits.length > 0) {
  console.error(
    `\n[check-no-banned-term] FAILED — banned term found in ${hits.length} location(s):\n`
  );
  for (const h of hits) console.error(`  ${h}`);
  console.error(
    '\nUse the approved phrasing: "effect-size-aware interval-null FDR gate ' +
      '(fdr < 0.05 at lfc = 0.25)".\n'
  );
  process.exit(1);
}

console.log(`[check-no-banned-term] OK — scanned ${files.length} file(s), no banned term.`);
