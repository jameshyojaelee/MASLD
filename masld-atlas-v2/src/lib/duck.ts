"use client";

/**
 * Generic DuckDB-WASM query layer for the portal.
 *
 * Thin, reusable helpers layered on the existing singleton in `./duckdb`
 * (`getDuckDB()`), so pages can register any Parquet in `public/data` and run
 * ad-hoc SQL against it without re-initialising a second database or worker.
 *
 * Transport is whole-file GET (fetch → arrayBuffer → `registerFileBuffer`),
 * NOT httpfs / HTTP range requests — a deliberate decision for the static
 * export + Hugging Face CDN deployment (see `data-base.ts`).
 *
 * Everything is lazy and client-only. Errors are surfaced to the caller
 * (rejected promises) rather than swallowed, so a page can render a fallback
 * instead of the whole route crashing.
 */

import type { AsyncDuckDB } from "@duckdb/duckdb-wasm";
import { getDuckDB } from "./duckdb";
import { dataUrl } from "./data-base";

// ---------------------------------------------------------------------------
// DB handle (reuse the existing singleton — never spin up a second worker)
// ---------------------------------------------------------------------------

/** Await/reuse the shared DuckDB instance created in `duckdb.ts`. */
export async function getDB(): Promise<AsyncDuckDB> {
  return getDuckDB();
}

// ---------------------------------------------------------------------------
// File registration (idempotent, de-duplicated across concurrent callers)
// ---------------------------------------------------------------------------

// Names already materialised as virtual files inside the WASM FS. The base
// singleton registers `atlas.parquet` during init, so seed it here to avoid a
// redundant (and potentially conflicting) re-registration.
const registeredFiles = new Set<string>(["atlas.parquet"]);
// In-flight registrations, keyed by virtual name, so parallel callers await one
// fetch rather than racing several downloads of the same buffer.
const pendingFiles = new Map<string, Promise<void>>();

/**
 * Register a Parquet file as a virtual file `name` inside DuckDB's FS.
 *
 * Idempotent: repeated calls for an already-registered `name` are no-ops, and
 * concurrent calls share a single fetch. After this resolves, the file can be
 * read in SQL via `read_parquet('<name>')`.
 */
export async function registerParquet(name: string, file: string): Promise<void> {
  if (registeredFiles.has(name)) return;
  const inflight = pendingFiles.get(name);
  if (inflight) return inflight;

  const task = (async () => {
    const db = await getDuckDB();
    const res = await fetch(dataUrl(file));
    if (!res.ok) {
      throw new Error(`registerParquet: failed to fetch ${file} (${res.status} ${res.statusText})`);
    }
    const buffer = new Uint8Array(await res.arrayBuffer());
    await db.registerFileBuffer(name, buffer);
    registeredFiles.add(name);
  })();

  // Track in-flight; clear the marker whether it resolves or rejects so a
  // failed fetch can be retried by a later call.
  pendingFiles.set(name, task);
  try {
    await task;
  } finally {
    pendingFiles.delete(name);
  }
}

// ---------------------------------------------------------------------------
// Row normalisation (Arrow → plain JS)
// ---------------------------------------------------------------------------

/** Convert an Arrow struct row to a plain object, coercing BigInt → Number. */
function normalizeRow(row: Record<string, unknown>): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const key in row) {
    const v = row[key];
    out[key] = typeof v === "bigint" ? Number(v) : v;
  }
  return out;
}

// ---------------------------------------------------------------------------
// Query execution
// ---------------------------------------------------------------------------

/**
 * Run SQL over a fresh connection and return plain JS rows.
 *
 * When `params` is provided the statement is prepared and bound (use `?`
 * placeholders) — the safe path for user-supplied values. Otherwise the SQL is
 * executed directly. The connection is always closed.
 */
export async function query<T = Record<string, unknown>>(
  sql: string,
  params?: unknown[]
): Promise<T[]> {
  const db = await getDuckDB();
  const conn = await db.connect();
  try {
    if (params && params.length > 0) {
      const stmt = await conn.prepare(sql);
      try {
        const result = await stmt.query(...params);
        return result.toArray().map((r: { toJSON: () => Record<string, unknown> }) =>
          normalizeRow(r.toJSON())
        ) as T[];
      } finally {
        await stmt.close();
      }
    }
    const result = await conn.query(sql);
    return result.toArray().map((r: { toJSON: () => Record<string, unknown> }) =>
      normalizeRow(r.toJSON())
    ) as T[];
  } finally {
    await conn.close();
  }
}

// ---------------------------------------------------------------------------
// SQL literal escaping
// ---------------------------------------------------------------------------

/**
 * Quote and escape a string as a SQL literal (doubles embedded single quotes).
 * Returns the value WITH surrounding quotes, e.g. `O'Neil` → `'O''Neil'`, so it
 * drops straight into an interpolated query. Prefer parameterised `query(sql,
 * params)` where possible; use this when building dynamic SQL by hand.
 */
export function sqlString(value: string): string {
  return `'${String(value).replace(/'/g, "''")}'`;
}

// ---------------------------------------------------------------------------
// Convenience: register a Parquet as a named table and query it
// ---------------------------------------------------------------------------

// Views created by `queryParquet`, so we only issue the DDL once per table.
const createdViews = new Set<string>();

/** Derive a safe SQL identifier (table name) from a file path's basename. */
function tableNameFromFile(file: string): string {
  const base = (file.split("/").pop() ?? file).replace(/\.[^.]+$/, "");
  const ident = base.replace(/[^A-Za-z0-9_]/g, "_");
  return /^[A-Za-z_]/.test(ident) ? ident : `t_${ident}`;
}

/**
 * Register `file` (basename → table name) and run `sql` against it.
 *
 * The Parquet is exposed as a view whose name is the file's basename without
 * extension, e.g. `atlas_core.parquet` → table `atlas_core`. Pass either a SQL
 * string that references that table, or a builder receiving the resolved table
 * name:
 *
 *   queryParquet("atlas_core.parquet",
 *     `SELECT * FROM atlas_core WHERE human_symbol = ${sqlString(sym)}`)
 *
 *   queryParquet("atlas_core.parquet",
 *     (t) => `SELECT * FROM ${t} WHERE human_symbol = ?`, [sym])   // parameterised
 */
export async function queryParquet<T = Record<string, unknown>>(
  file: string,
  sql: string | ((table: string) => string),
  params?: unknown[]
): Promise<T[]> {
  const table = tableNameFromFile(file);
  const fileName = file.split("/").pop() ?? file;

  await registerParquet(fileName, file);

  if (!createdViews.has(table)) {
    const db = await getDuckDB();
    const conn = await db.connect();
    try {
      await conn.query(
        `CREATE VIEW IF NOT EXISTS "${table}" AS SELECT * FROM read_parquet('${fileName}')`
      );
      createdViews.add(table);
    } finally {
      await conn.close();
    }
  }

  const text = typeof sql === "function" ? sql(table) : sql;
  return query<T>(text, params);
}
