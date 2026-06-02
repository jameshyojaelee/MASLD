import * as duckdb from "@duckdb/duckdb-wasm";

let db: duckdb.AsyncDuckDB | null = null;
let initPromise: Promise<duckdb.AsyncDuckDB> | null = null;

export async function getDuckDB(): Promise<duckdb.AsyncDuckDB> {
  if (db) return db;
  if (initPromise) return initPromise;

  initPromise = (async () => {
    const JSDELIVR_BUNDLES = duckdb.getJsDelivrBundles();
    const bundle = await duckdb.selectBundle(JSDELIVR_BUNDLES);

    const worker_url = URL.createObjectURL(
      new Blob([`importScripts("${bundle.mainWorker!}");`], {
        type: "text/javascript",
      })
    );

    const worker = new Worker(worker_url);
    const logger = new duckdb.ConsoleLogger();
    db = new duckdb.AsyncDuckDB(logger, worker);
    await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
    URL.revokeObjectURL(worker_url);

    const response = await fetch("/data/atlas.parquet");
    const buffer = await response.arrayBuffer();
    await db.registerFileBuffer("atlas.parquet", new Uint8Array(buffer));

    const conn = await db.connect();
    await conn.query(
      "CREATE VIEW atlas AS SELECT * FROM read_parquet('atlas.parquet')"
    );
    await conn.close();

    return db;
  })();

  return initPromise;
}

export async function queryAtlas<T = Record<string, unknown>>(
  sql: string
): Promise<T[]> {
  const instance = await getDuckDB();
  const conn = await instance.connect();
  try {
    const result = await conn.query(sql);
    return result.toArray().map((row: any) => row.toJSON()) as T[];
  } finally {
    await conn.close();
  }
}

export async function getAtlasRowCount(): Promise<number> {
  const rows = await queryAtlas<{ cnt: number }>(
    "SELECT COUNT(*) as cnt FROM atlas"
  );
  return rows[0]?.cnt ?? 0;
}
