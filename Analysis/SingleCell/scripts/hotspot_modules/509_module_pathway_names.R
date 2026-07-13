#!/usr/bin/env Rscript
# ============================================================================
# 509_module_pathway_names.R
#
# Functional NAME for every Hotspot module via over-representation analysis
# (hypergeometric) of each module's member genes against NAMED MSigDB
# collections (Hallmark + Reactome + GO-BP), using that cell type's Hotspot
# autocorrelation gene set as the BACKGROUND UNIVERSE (not the whole genome —
# the ~2000 autocorr genes are already a biased HVG-like set, so a genome-wide
# background would just recover "highly-expressed gene" pathways).
#
# A module is a discrete weighted gene SET, so ORA (not phenotype-ranked GSEA)
# is the appropriate test. Names prefer Hallmark (50 broad, clean sets); the
# single most-significant named pathway across all collections is also kept.
#
# Outputs (authoritative, survive a 508 rebuild):
#   module_pathway_enrichment.tsv   full ORA table (all sig pathways per module)
#   module_names.tsv                cell_type, module, module_name, hallmark, top_pathway
# Convenience: merges module_name / module_hallmark / module_top_pathway into
#   all_modules.tsv (NOTE: 508 rebuilds all_modules.tsv from scratch — re-run
#   this script after any 508 to restore the columns).
#
# 95d-style hand-curation: edit CURATED below (key "celltype__module") to
# override the auto name; re-run to apply.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(msigdbr) })

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HS <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
AM <- file.path(HS, "all_modules.tsv")
CTS <- c("hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "tcells")

MIN_SET   <- 10   # min pathway genes (after universe intersect)
MIN_MODULE <- 5   # min module genes (after universe intersect)
MIN_OVERLAP <- 2

# Hand-curated names for ALL modules (95d pattern); key = "<cell_type>__<module>".
# Approved 2026-07-01 ("curate all"): gene-grounded names disambiguating the
# Hallmark-name collisions (e.g. the 5 hepatocyte "Coagulation" modules are
# distinct — apolipoprotein / acute-phase / albumin-secretory / coag-factor /
# complement). Reinterpreted from the auto Hallmark where the top-weighted genes
# disagree (e.g. hep__24 TNFa->NRF2 antioxidant; hep__2 KRAS->sinusoidal endo).
CURATED <- list(
  # ---- Hepatocytes ----
  "hepatocytes__1"  = "Leukocyte / immune (DOCK2)",
  "hepatocytes__2"  = "Sinusoidal endothelial (STAB2)",
  "hepatocytes__3"  = "Oxidative phosphorylation",
  "hepatocytes__4"  = "Cell adhesion (CNTNAP2)",
  "hepatocytes__5"  = "Apolipoprotein / plasma (APOA2)",
  "hepatocytes__6"  = "Hypoxia (NDRG1)",
  "hepatocytes__7"  = "Acute-phase / complement (C3)",
  "hepatocytes__8"  = "Stromal ECM (IGFBP7)",
  "hepatocytes__9"  = "Noncoding RNA",
  "hepatocytes__10" = "Lipid droplet (PLIN5)",
  "hepatocytes__11" = "Mitochondrial ETC (NDUF/COX)",
  "hepatocytes__12" = "Secretory plasma protein (ALB)",
  "hepatocytes__13" = "Nucleoside catabolism",
  "hepatocytes__14" = "Amino-acid metabolism (GOT1)",
  "hepatocytes__15" = "Coagulation factor (F12)",
  "hepatocytes__16" = "PPARa lipid metabolism",
  "hepatocytes__17" = "Hepatocyte stress (SLC38A1)",
  "hepatocytes__18" = "Amino-acid degradation (TAT)",
  "hepatocytes__19" = "Fatty-acid / peroxisomal metab",
  "hepatocytes__20" = "Ductular injury (BICC1)",
  "hepatocytes__21" = "Xenobiotic / drug metab (CYP)",
  "hepatocytes__22" = "Ion / membrane potential",
  "hepatocytes__23" = "IL6-STAT3 acute phase (CRP)",
  "hepatocytes__24" = "NRF2 antioxidant (TXNRD1)",
  "hepatocytes__25" = "UPR / ERN1 stress",
  "hepatocytes__26" = "AP-1 stress (ATF3/GDF15)",
  "hepatocytes__27" = "Complement (CFH)",
  "hepatocytes__28" = "Lipid / autophagy (VMP1)",
  "hepatocytes__29" = "Cell polarity (PARD3)",
  "hepatocytes__30" = "Immediate-early AP-1 (FOS)",
  # ---- Cholangiocytes ----
  "cholangiocytes__1"  = "Plasma lipoprotein (APOE)",
  "cholangiocytes__2"  = "Oxidative phosphorylation",
  "cholangiocytes__3"  = "Macrophage / immune (CD163)",
  "cholangiocytes__4"  = "Cholangiocyte secretory (TFF3)",
  "cholangiocytes__5"  = "Mitochondrial / metabolic (ETFB)",
  "cholangiocytes__6"  = "Noncoding RNA",
  "cholangiocytes__7"  = "Uncharacterized (DACH1)",
  "cholangiocytes__8"  = "Nutrient sensing (CASR)",
  "cholangiocytes__9"  = "Junction / autophagy (TJP1)",
  "cholangiocytes__10" = "Complement / plasma (CFH)",
  "cholangiocytes__11" = "Ductal morphogenesis (ANXA4)",
  "cholangiocytes__12" = "Ciliary / apical (PKHD1)",
  "cholangiocytes__13" = "Bile-acid transport (ABCB11)",
  "cholangiocytes__14" = "Cell adhesion / actin",
  "cholangiocytes__15" = "Hypoxia (NDRG1)",
  "cholangiocytes__16" = "Mixed metabolic (CYP2C9)",
  "cholangiocytes__17" = "cAMP / PDE signaling",
  "cholangiocytes__18" = "CFTR / ion transport",
  "cholangiocytes__19" = "Cytoskeleton / TGFB (MAP4)",
  "cholangiocytes__20" = "Glucocorticoid / stress (FKBP5)",
  "cholangiocytes__21" = "Inflammatory / oxidative (SOD2)",
  "cholangiocytes__22" = "Progenitor / Notch (MSI2)",
  "cholangiocytes__23" = "Mesenchymal / EMT (VIM)",
  "cholangiocytes__24" = "Cytoskeleton / FGFR2",
  "cholangiocytes__25" = "FXR / bile-acid (NR1H4)",
  "cholangiocytes__26" = "Adhesion / polarity (ALCAM)",
  "cholangiocytes__27" = "NF-kB immune activation",
  "cholangiocytes__28" = "Rho-GTPase / adhesion",
  # ---- Fibroblasts ----
  "fibroblasts__1"  = "Plasma protein (APOA1)",
  "fibroblasts__2"  = "Oxidative phosphorylation",
  "fibroblasts__3"  = "Cell adhesion / junction",
  "fibroblasts__4"  = "Portal fibroblast (NR1H4)",
  "fibroblasts__5"  = "Myofibroblast (EBF1/SLIT3)",
  "fibroblasts__6"  = "Biliary-like (PKHD1/ANXA4)",
  "fibroblasts__7"  = "Contractile myofibroblast",
  "fibroblasts__8"  = "Lipogenic hep-like (MLXIPL)",
  "fibroblasts__9"  = "Hypoxia (NDRG1)",
  "fibroblasts__10" = "Hepatocyte xenobiotic (CYP2C)",
  "fibroblasts__11" = "cAMP / PKA signaling",
  "fibroblasts__12" = "S100 / mesenchymal (VIM)",
  "fibroblasts__13" = "Myeloid / immune (DOCK2)",
  "fibroblasts__14" = "Activated stellate (PDGFRA)",
  "fibroblasts__15" = "Hypoxia / VEGF (HILPDA)",
  "fibroblasts__16" = "Actomyosin / MHC-I",
  "fibroblasts__17" = "Elastic-fiber ECM (FBLN5)",
  "fibroblasts__18" = "Cell motility / filopodium",
  "fibroblasts__19" = "ECM proteoglycan (BGN)",
  "fibroblasts__20" = "Neuronal adhesion (GRID2)",
  "fibroblasts__21" = "Cytoskeleton / ECM (NPNT)",
  "fibroblasts__22" = "Adhesion / signaling (EPHA3)",
  "fibroblasts__23" = "Insulin response (FOXO1)",
  "fibroblasts__24" = "Activated fibroblast (MMP19)",
  "fibroblasts__25" = "Notch / RXR (NOTCH1)",
  "fibroblasts__26" = "Versican / GAG synthesis",
  "fibroblasts__27" = "TGFB / collagen (ADAMTS2)",
  "fibroblasts__28" = "ECM adhesion (TNC)",
  "fibroblasts__29" = "ECM organization",
  # ---- Macrophages ----
  "macrophages__1"  = "Myeloid housekeeping (TYROBP)",
  "macrophages__2"  = "Plasma protein (APOA1)",
  "macrophages__3"  = "Glycolysis / hypoxia (HK2)",
  "macrophages__4"  = "Adhesion (ITGA9/TCF7L2)",
  "macrophages__5"  = "Acute-phase / MT (SAA1)",
  "macrophages__6"  = "Vesicle / signaling (ITSN1)",
  "macrophages__7"  = "Glucocorticoid resp (FKBP5)",
  "macrophages__8"  = "Myeloid differentiation (DOCK2)",
  "macrophages__9"  = "Inflammatory myeloid (TNFAIP2)",
  "macrophages__10" = "Chromatin / epigenetic (DNMT3A)",
  "macrophages__11" = "Vesicle transport (VPS13B)",
  "macrophages__12" = "Resident macrophage (LIPA)",
  "macrophages__13" = "Actin cytoskeleton (DOCK8)",
  "macrophages__14" = "Adhesion / Rho (ARHGAP22)",
  "macrophages__15" = "Endocytosis / TLR2",
  "macrophages__16" = "Hep-like metabolic (CYP3A5)",
  "macrophages__17" = "Membrane dynamics (DNM2)",
  "macrophages__18" = "Lipid-associated Mac (PLIN5)",
  "macrophages__19" = "Inflammatory activation (CD83)",
  # ---- T cells ----
  "tcells__1"  = "Plasma protein (APOC1)",
  "tcells__2"  = "T-cell signaling (PTPRC)",
  "tcells__3"  = "Cytotoxic / MHC-I (HLA-B)",
  "tcells__4"  = "Translation / ribosome (EEF1A1)",
  "tcells__5"  = "T-cell activation (PIK3CD)",
  "tcells__6"  = "Rho-GTPase (DOCK11)",
  "tcells__7"  = "Chromatin remodeling (KDM2A)",
  "tcells__8"  = "Memory / MAIT (IL7R)",
  "tcells__9"  = "NK / cytotoxic (KLRD1)",
  "tcells__10" = "Translation / MYC (EEF1)",
  "tcells__11" = "Immediate-early (ZFP36/JUNB)"
)

# ---- named collections only (skip cNMF/NMF/SCENIC which yield 'P#') ----
message("[msigdbr] loading Hallmark + Reactome + GO-BP ...")
grab <- function(coll, sub = NULL) {
  a <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
       else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  as.data.table(a)[, .(gs_name, gene_symbol)]
}
msig <- rbindlist(list(
  grab("H")[,                     coll := "HALLMARK"],
  grab("C2", "CP:REACTOME")[,     coll := "REACTOME"],
  grab("C5", "GO:BP")[,           coll := "GOBP"]
))
message(sprintf("[msigdbr] %d gene-set memberships across 3 collections",
                nrow(msig)))

clean_name <- function(gs) {
  x <- sub("^(HALLMARK|REACTOME|GOBP|GOMF|GOCC|KEGG|WP)_", "", gs)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

ora_ct <- function(ct) {
  mgf <- file.path(HS, ct, "module_genes.tsv")
  acf <- file.path(HS, ct, "autocorr.tsv")
  if (!file.exists(mgf) || !file.exists(acf)) {
    message(sprintf("[skip] %s (missing module_genes/autocorr)", ct)); return(NULL)
  }
  mg   <- fread(mgf)
  univ <- unique(fread(acf)$gene)
  N    <- length(univ)
  # restrict pathways to the universe, drop tiny sets
  msu   <- msig[gene_symbol %in% univ]
  paths <- split(msu$gene_symbol, paste(msu$coll, msu$gs_name, sep = "|"))
  paths <- paths[lengths(paths) >= MIN_SET]
  message(sprintf("[%s] universe=%d genes, %d pathways after universe filter",
                  ct, N, length(paths)))
  pk <- lengths(paths)
  res <- rbindlist(lapply(sort(unique(mg$module)), function(m) {
    if (m == 0) return(NULL)
    genes <- intersect(mg[module == m, gene], univ)
    mlen  <- length(genes)
    if (mlen < MIN_MODULE) return(NULL)
    a <- vapply(paths, function(pg) length(intersect(genes, pg)), integer(1))
    keep <- a >= MIN_OVERLAP
    if (!any(keep)) return(NULL)
    k <- pk[keep]; av <- a[keep]
    p <- phyper(av - 1, k, N - k, mlen, lower.tail = FALSE)
    nm <- names(paths)[keep]
    data.table(cell_type = ct, module = m,
               coll = sub("\\|.*", "", nm), gs = sub(".*\\|", "", nm),
               overlap = av, path_size = k, mod_size = mlen, p = p)
  }))
  if (!nrow(res)) return(NULL)
  res[, padj := p.adjust(p, "BH"), by = module]      # per-module BH
  res
}

all_ora <- rbindlist(lapply(CTS, ora_ct))
setorder(all_ora, cell_type, module, padj, p)
fwrite(all_ora, file.path(HS, "module_pathway_enrichment.tsv"), sep = "\t")

# ---- per-module top annotation ----
pick <- function(d) {
  hall <- d[coll == "HALLMARK" & padj < 0.05][order(padj, p)][1]
  top  <- d[padj < 0.05][order(padj, p)][1]
  data.table(
    module_hallmark = if (nrow(hall) && !is.na(hall$gs))
      sprintf("%s (q=%.1e)", clean_name(hall$gs), hall$padj) else NA_character_,
    module_top_pathway = if (nrow(top) && !is.na(top$gs))
      sprintf("%s: %s (q=%.1e)", top$coll, clean_name(top$gs), top$padj) else NA_character_,
    module_name = if (nrow(hall) && !is.na(hall$gs)) clean_name(hall$gs)
                  else if (nrow(top) && !is.na(top$gs)) clean_name(top$gs)
                  else "Unresolved")
}
ann <- all_ora[, pick(.SD), by = .(cell_type, module)]
# apply hand-curated overrides
ann[, key := paste(cell_type, module, sep = "__")]
for (k in names(CURATED)) ann[key == k, module_name := CURATED[[k]]]
ann[, key := NULL]
fwrite(ann, file.path(HS, "module_names.tsv"), sep = "\t")

# ---- merge into all_modules.tsv (convenience; clobbered by 508 re-run) ----
am <- fread(AM); am[, module := as.integer(module)]
drop <- intersect(c("module_hallmark", "module_top_pathway", "module_name"), names(am))
if (length(drop)) am[, (drop) := NULL]
am <- merge(am, ann, by = c("cell_type", "module"), all.x = TRUE)
am[is.na(module_name), module_name := "Unresolved"]
fwrite(am, AM, sep = "\t")

cat(sprintf("\nAnnotated %d modules across %d cell types.\n",
            nrow(ann), uniqueN(ann$cell_type)))
cat("Resolved (named) vs Unresolved:\n")
print(ann[, .N, by = .(named = module_name != "Unresolved")])
cat("\nPreview (hep + a few per CT):\n")
print(ann[order(cell_type, module), .(cell_type, module, module_name)])
