#!/usr/bin/env Rscript
# AN5 task 3: lead-LD rule (r2 >= 0.8) at the caQTL-eQTL step of the H53 three-way link.
# Reads H53 outputs and inputs read-only; writes only under OUT. Deterministic (no random draws; seed set anyway).
# eQTL lead   = top-alpha variant of the colocalized eQTL SuSiE effect (idx2) in results/eqtl_susie_polyfun.
# caQTL lead  = variant with the smallest caQTL p among positions shared between the eQTL row and the
#               caQTL peak row, keyed exactly as h53_job2.R:242-249 (the same rule h53_common.R:202-205
#               uses for the GWAS-caQTL step).
# LD          = PolyFun UK Biobank EUR panel (data/ld_ref/polyfun_eur), the panel the eQTL SuSiE fits used
#               (src/08_run_eqtl_susie.R, LD_PANEL=polyfun). Blocks are assembled block-diagonally in that
#               fine-mapping, so two variants in different approx-LD blocks have r = 0 in the panel as used.
set.seed(20260923)
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z/scripts/h53_common.R")
OUT <- Sys.getenv("AN5_OUT"); stopifnot(nzchar(OUT))
H53 <- RUN
LDB <- file.path(FM, "data/ld_ref/polyfun_eur")
EQ  <- file.path(PROJ, "data/broadaway_eqtl")

TW  <- fread(file.path(H53, "results/three_way.tsv"))
G   <- fread(file.path(H53, "results/catalog_gene_caqtl_annotation.tsv"))
SL  <- fread(file.path(H53, "results/study_leads_job1.tsv"), colClasses = list(character = "chr"))
SL[, lead_key := paste(gwas, hit1, sep = "|")]
a   <- load_arms()
Lg  <- unique(a$L[, .(lead_key = paste(gwas, hit1, sep = "|"), ensembl, gene, idx2)])

M1_printed <- G[caqtl_three_way_eqtl == TRUE, unique(ensembl)]
stamp("M1 printed:", length(M1_printed))
Q <- TW[!is.na(tw_PP.H4) & tw_PP.H4 > PP_ENTRY]
Q <- merge(Q, unique(Lg[, .(lead_key, ensembl)]), by = c("lead_key", "ensembl"))   # liver-arm gene rows only
Q <- Q[ensembl %in% M1_printed]
stopifnot(uniqueN(Q$ensembl) == length(M1_printed))
Q[, chr := sub(":.*", "", sub("^[^|]*\\|", "", lead_key))]
Q <- merge(Q, SL[, .(lead_key, refit_ensembl, refit_gene, gwas_lead_snp)], by = "lead_key", all.x = TRUE)
stamp("qualifying three-way rows:", nrow(Q), "genes:", uniqueN(Q$ensembl))

blocks <- fread(file.path(LDB, "approx_LD_blocks.txt"))
block_of <- function(ch, pos) {
  b <- blocks[chr == as.integer(ch) & start <= pos & stop > pos]
  if (!nrow(b)) return(NA_character_)
  file.path(LDB, paste0("chr", ch), paste0(b$start[1], ".", b$stop[1]), paste0(b$start[1], ".", b$stop[1]))
}
.bim <- new.env()
get_bim <- function(f) {
  if (is.null(.bim[[f]])) { b <- fread(paste0(f, ".bim"), header = FALSE); setnames(b, c("chr", "rsid", "cm", "pos", "a1", "a2")); b[, line := .I]; .bim[[f]] <- b }
  .bim[[f]]
}
bim_line <- function(f, p0, A1, A2) {
  b <- get_bim(f); h <- b[pos == p0]
  if (!nrow(h)) return(NA_integer_)
  m <- h[(a1 == A1 & a2 == A2) | (a1 == A2 & a2 == A1)]
  if (nrow(m)) m$line[1] else NA_integer_
}

rows <- list()
for (ch in sort(unique(Q$chr))) {
  qc <- Q[chr == ch]
  stamp("chr", ch, ":", nrow(qc), "rows")
  f <- list.files(CQDIR, pattern = paste0("^caQTL_hg38_chr", ch, "_.*txt.gz$"), full.names = TRUE)
  q <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", showProgress = FALSE); setkey(q, peak)
  mq <- fread(file.path(EQ, paste0("chr", ch, "_marginal_summary_results.tsv")),
              select = c("POS", "NEA", "EA", "PVAL", "ENSG"), colClasses = list(character = "PVAL"))
  mq <- mq[ENSG %in% qc$ensembl]
  mq[, log10p := { x <- tolower(PVAL); m <- suppressWarnings(as.numeric(sub("e.*", "", x))); e <- suppressWarnings(as.numeric(sub(".*e", "", x)))
                   ifelse(grepl("e", x), log10(m) + e, log10(as.numeric(x))) }]
  for (i in seq_len(nrow(qc))) {
    r <- qc[i]
    e <- readRDS(file.path(EQTL_SUSIE_DIR, paste0("chr", ch), paste0(r$ensembl, "_susie.rds")))
    ids <- colnames(e$lbf_variable)
    elead <- ids[which.max(e$alpha[r$idx2, ])]
    p38 <- lift19to38(ch, as.integer(sub("^[0-9]+:", "", ids)))
    k <- ifelse(is.na(p38), paste0("e|", ids), paste0("p", p38)); d <- duplicated(k); if (any(d)) k[d] <- paste0("edup|", ids[d])
    qp <- caqtl_prepare(copy(q[.(r$peak)])); qp <- qp[order(pvalue)]
    kpos <- paste0("p", qp$pos38); dq <- duplicated(kpos); kpos[dq] <- paste0("qdup|", qp$variant[dq])
    sh <- intersect(kpos, k)
    cl_key <- sh[1]  # qp is p-ordered, intersect keeps order of kpos
    cl_var <- qp$variant[match(cl_key, kpos)]; cl_p <- qp$pvalue[match(cl_key, kpos)]
    clead <- ids[match(cl_key, k)]
    # marginal eQTL support inside the idx2 credible set
    csn <- paste0("L", r$idx2)
    cs <- e$sets$cs[[csn]]
    cs_pos <- if (is.null(cs)) integer(0) else as.integer(sub("^[0-9]+:", "", ids[cs]))
    mg <- mq[ENSG == r$ensembl]
    cs_min_log10p <- if (length(cs_pos)) suppressWarnings(min(mg[POS %in% cs_pos, log10p])) else NA_real_
    # LD in the PolyFun EUR panel
    pe <- as.integer(sub("^[0-9]+:", "", elead)); pc <- as.integer(sub("^[0-9]+:", "", clead))
    ae <- mg[POS == pe][1]; ac <- mg[POS == pc][1]
    fe <- block_of(ch, pe); fc <- block_of(ch, pc)
    le <- if (!is.na(fe)) bim_line(fe, pe, ae$NEA, ae$EA) else NA_integer_
    lc <- if (!is.na(fc)) bim_line(fc, pc, ac$NEA, ac$EA) else NA_integer_
    rows[[length(rows) + 1L]] <- data.table(r[, .(lead_key, peak, gene, ensembl, idx2, tw_PP.H4, r2_eqtl_pp4, refit_ensembl, refit_gene)],
      eqtl_lead19 = elead, eqtl_lead_alpha = max(e$alpha[r$idx2, ]), caqtl_lead_variant38 = cl_var, caqtl_lead_p = cl_p,
      caqtl_lead19 = clead, same_variant = identical(elead, clead), n_shared = length(sh),
      cs_size = length(cs), cs_min_marginal_log10p = cs_min_log10p,
      block_e = fe, block_c = fc, line_e = le, line_c = lc)
  }
  rm(q, mq); gc(FALSE)
}
R <- rbindlist(rows)
# extract r with one awk pass per block (rows needed: line_e, field line_c)
R[, r_ld := NA_real_]
R[same_variant == TRUE, r_ld := 1]
R[!same_variant & !is.na(block_e) & !is.na(block_c) & block_e != block_c, r_ld := 0]
need <- R[!same_variant & block_e == block_c & !is.na(line_e) & !is.na(line_c)]
for (bf in unique(need$block_e)) {
  nb <- need[block_e == bf]
  tmp <- tempfile(); fwrite(unique(nb[, .(line_e, line_c)]), tmp, sep = "\t", col.names = FALSE)
  cmd <- sprintf("awk -v F=%s 'BEGIN{while((getline l < F)>0){split(l,a,\"\\t\"); w[a[1]]=w[a[1]] \" \" a[2]; if(a[1]+0>mx) mx=a[1]+0}} (NR in w){n=split(w[NR],c,\" \"); for(k=1;k<=n;k++) print NR\"\\t\"c[k]\"\\t\"$(c[k])} NR>=mx{exit}' %s",
                 shQuote(tmp), shQuote(paste0(bf, ".ld")))
  v <- fread(cmd = cmd, header = FALSE, col.names = c("line_e", "line_c", "r"))
  for (j in seq_len(nrow(v))) R[block_e == bf & line_e == v$line_e[j] & line_c == v$line_c[j], r_ld := v$r[j]]
  unlink(tmp)
}
R[, r2_ld := r_ld^2]
R[, ld_status := fifelse(same_variant, "same variant",
                 fifelse(is.na(block_e) | is.na(block_c) | (block_e == block_c & (is.na(line_e) | is.na(line_c))), "not in panel",
                 fifelse(block_e != block_c, "different LD blocks (r = 0 in panel as used)", "same block")))]
fwrite(R, file.path(OUT, "threeway_rows_lead_ld.tsv"), sep = "\t")

gene <- R[, .(gene = gene[1], n_rows = .N, best_tw_pp4 = max(tw_PP.H4),
              any_r2_ge_0.8 = any(r2_ld >= 0.8, na.rm = TRUE),
              any_r2_ge_0.5 = any(r2_ld >= 0.5, na.rm = TRUE),
              max_r2 = suppressWarnings(max(r2_ld, na.rm = TRUE)),
              n_rows_unresolved = sum(is.na(r2_ld)),
              own_refit_row = any(refit_ensembl == ensembl),
              best_cs_min_marginal_log10p = suppressWarnings(min(cs_min_marginal_log10p, na.rm = TRUE))), by = ensembl]
fwrite(gene, file.path(OUT, "threeway_gene_lead_ld.tsv"), sep = "\t")
summ <- list(M1_printed = length(M1_printed), qualifying_rows = nrow(R),
             M1_under_r2_ge_0.8 = gene[any_r2_ge_0.8 == TRUE, .N],
             M1_under_r2_ge_0.5 = gene[any_r2_ge_0.5 == TRUE, .N],
             genes_with_unresolved_rows_only = gene[n_rows_unresolved == n_rows, .N],
             rows_by_ld_status = as.list(table(R$ld_status)),
             genes_no_own_refit_row = gene[own_refit_row == FALSE, sort(gene)],
             genes_cs_marginal_p_gt_1e4 = gene[best_cs_min_marginal_log10p > -4, sort(gene)],
             genes_failing_r2_0.8 = gene[any_r2_ge_0.8 == FALSE, sort(gene)])
jsonlite::write_json(summ, file.path(OUT, "threeway_summary.json"), auto_unbox = TRUE, pretty = TRUE)
print(summ)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
stamp("done")
