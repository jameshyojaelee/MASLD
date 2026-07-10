import puppeteer from "puppeteer";
const BASE = "https://jameshyojaelee-masld-atlas.static.hf.space";
const DIR = "/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/144cfaac-c7b7-41db-ba80-2612be2de99e/scratchpad";
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

console.log("waiting 25s for propagation…");
await sleep(25000);
const browser = await puppeteer.launch({ headless: "new", args: ["--no-sandbox", "--disable-setuid-sandbox"] });

async function prep(page) {
  await page.evaluate(async () => {
    for (let y = 0; y <= document.body.scrollHeight; y += 300) { window.scrollTo(0, y); await new Promise((r) => setTimeout(r, 90)); }
    window.scrollTo(0, 0);
  });
  await sleep(2200); // let UMAP bloom settle + a couple drift frames
}

async function assert(page) {
  return page.evaluate(() => {
    const h2s = [...document.querySelectorAll("h2")].map((h) => h.textContent.trim());
    return {
      dark: document.documentElement.classList.contains("dark"),
      kpi_cards_gone: !document.body.textContent.includes("COLOC EFFECTORS") && !document.body.textContent.includes("CONVERGENT"),
      at_a_glance_gone: !h2s.includes("At a glance"),
      has_bento: h2s.includes("Explore the atlas"),
      has_featured: h2s.includes("Featured genes"),
      umap_canvas: !!document.querySelector('canvas[aria-label="Single-cell UMAP thumbnail"]'),
      aurora_present: !!document.querySelector('[class*="aurora"], .aurora-blob') || document.querySelectorAll("section canvas, section div[aria-hidden]").length > 0,
      data_error: /Entry not found|is not valid JSON/.test(document.body.textContent || ""),
    };
  });
}

// Dark (default)
const d = await browser.newPage();
await d.setViewport({ width: 1440, height: 1000, deviceScaleFactor: 2 });
await d.goto(`${BASE}/?cb=${Date.now()}#/`, { waitUntil: "networkidle0", timeout: 60000 });
await sleep(3500);
await prep(d);
console.log("DARK", JSON.stringify(await assert(d)));
await d.screenshot({ path: `${DIR}/home_r3_dark.png`, fullPage: true });
// hero-only crop for aurora/constellation detail
await d.screenshot({ path: `${DIR}/home_r3_hero.png`, clip: { x: 210, y: 0, width: 1230, height: 640 } });

// Light
await d.evaluate(() => localStorage.setItem("masld-atlas-theme", "light"));
await d.reload({ waitUntil: "networkidle0", timeout: 60000 });
await sleep(3500);
await prep(d);
console.log("LIGHT", JSON.stringify(await assert(d)));
await d.screenshot({ path: `${DIR}/home_r3_light.png`, fullPage: true });

await browser.close();
console.log("shots saved");
