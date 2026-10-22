          ## ==== AN5 diagnostic block (REVIEW2 triage 2026-10-21; inserted into a COPY of src/06_susie_coloc.R) ====
          ## Records the GWAS-side SuSiE fit, a kriging_rss z-vs-LD consistency check, and three refits:
          ## (b) LD diagonal 1e-6 instead of 1e-3, (c) residual variance estimated, (d) kriging outliers removed.
          ## Each fit is colocalized against the same full eQTL fit with the same eligibility rule.
          tryCatch({
            an5_dir <- Sys.getenv("AN5_DIAG_DIR")
            an5_tag <- paste(gene_symbol, gwas_name, sep = "__")
            an5_fit <- function(Rm, erv = FALSE, keep = rep(TRUE, length(gwas_z))) {
              zz <- gwas_z[keep]; RR <- Rm[keep, keep, drop = FALSE]
              args <- list(z = zz, R = RR, n = gwas_n_eff, L = SUSIE_L,
                           estimate_residual_variance = erv, check_R = FALSE,
                           min_abs_corr = 0.5, max_iter = 500)
              if (!erv) args$residual_variance <- 1
              set.seed(42)
              tryCatch(do.call(susie_rss, args), error = function(e) {
                cat("  AN5 fit error:", conditionMessage(e), "\n"); NULL })
            }
            an5_summ <- function(fit, label, zz) {
              if (is.null(fit)) return(data.table(fit = label, fit_ok = FALSE))
              cs <- fit$sets$cs; ncs <- length(cs)
              sz <- if (ncs) vapply(cs, length, 1L) else integer(0)
              pur <- if (ncs) fit$sets$purity$min.abs.corr else numeric(0)
              ld <- if (ncs) vapply(cs, function(ix) names(zz)[ix][which.max(fit$pip[ix])], "") else character(0)
              lz <- if (ncs) vapply(cs, function(ix) zz[ix][which.max(fit$pip[ix])], 0) else numeric(0)
              data.table(fit = label, fit_ok = TRUE, converged = isTRUE(fit$converged),
                         niter = fit$niter, n_cs = ncs,
                         cs_sizes = paste(sz, collapse = ";"),
                         cs_min_purity = paste(signif(pur, 3), collapse = ";"),
                         cs_lead = paste(ld, collapse = ";"),
                         cs_lead_z = paste(round(lz, 2), collapse = ";"),
                         sigma2 = as.numeric(fit$sigma2))
            }
            an5_coloc <- function(fit, label) {
              empty <- data.table(fit = label, best_pp4 = NA_real_, best_pp3 = NA_real_,
                                  n_pairs = NA_integer_, n_eligible = NA_integer_,
                                  n_pp4_gt_0.5 = NA_integer_, best_hit1 = NA_character_, best_hit2 = NA_character_)
              if (is.null(fit) || !isTRUE(fit$converged) || length(fit$sets$cs) == 0) return(list(s = empty, p = NULL))
              set.seed(42)
              res <- tryCatch(coloc.susie(fit, s_eqtl), error = function(e) NULL)
              if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0 ||
                  !"PP.H4.abf" %in% names(res$summary)) return(list(s = empty, p = NULL))
              pt <- coloc_susie_pair_table(fit, s_eqtl, res)
              kp <- pt$pairs[eligible %in% TRUE & is.finite(PP.H4)]
              if (nrow(kp) == 0) { empty$n_pairs <- nrow(pt$pairs); empty$n_eligible <- 0L; return(list(s = empty, p = pt$pairs)) }
              b <- kp[which.max(PP.H4)]
              list(s = data.table(fit = label, best_pp4 = b$PP.H4, best_pp3 = b$PP.H3,
                                  n_pairs = nrow(pt$pairs), n_eligible = nrow(kp),
                                  n_pp4_gt_0.5 = sum(kp$PP.H4 > 0.5),
                                  best_hit1 = b$hit1, best_hit2 = b$hit2),
                   p = cbind(fit = label, pt$pairs))
            }
            # (a) z-vs-LD consistency: susieR kriging_rss on the 1e-3 regularised panel
            kr <- tryCatch(susieR::kriging_rss(z = gwas_z, R = R_reg, n = gwas_n_eff),
                           error = function(e) { cat("  AN5 kriging error:", conditionMessage(e), "\n"); NULL })
            flag <- rep(FALSE, length(gwas_z))
            kr_row <- data.table(n_snps = length(gwas_z), kriging_ok = !is.null(kr))
            if (!is.null(kr)) {
              cd <- as.data.table(kr$conditional_dist)
              flag <- (cd$logLR > 2 & abs(cd$z) > 2)
              flag[is.na(flag)] <- FALSE
              top <- which.max(abs(gwas_z))
              kr_row <- data.table(n_snps = length(gwas_z), kriging_ok = TRUE, n_dentist5.45 = sum(abs(cd$z_std_diff) > 5.45, na.rm = TRUE),
                                   n_flag_logLR2_absz2 = sum(flag),
                                   max_abs_z_std_diff = max(abs(cd$z_std_diff), na.rm = TRUE),
                                   n_abs_z_std_diff_gt3 = sum(abs(cd$z_std_diff) > 3, na.rm = TRUE),
                                   lead_snp = names(gwas_z)[top], lead_z = gwas_z[top],
                                   lead_condmean = cd$condmean[top], lead_z_std_diff = cd$z_std_diff[top],
                                   n_absz_gt5 = sum(abs(gwas_z) > 5),
                                   n_flag_among_absz_gt5 = sum(flag & abs(gwas_z) > 5))
              }
            # DENTIST-equivalent outlier rule: |z - E[z | LD]| / sd > 5.45 (chi-square 1 df, p < 5e-8)
            dflag <- rep(FALSE, length(gwas_z))
            if (!is.null(kr)) { dflag <- abs(cd$z_std_diff) > 5.45; dflag[is.na(dflag)] <- FALSE }
            fwrite(data.table(snp = names(gwas_z)[dflag]), file.path(an5_dir, paste0(an5_tag, "__dentist_flagged.tsv")), sep = "\t")
            an5_fitL <- function(Rm, Lval) {
              set.seed(42)
              tryCatch(susie_rss(z = gwas_z, R = Rm, n = gwas_n_eff, L = Lval, estimate_residual_variance = FALSE,
                                 residual_variance = 1, check_R = FALSE, min_abs_corr = 0.5, max_iter = 500),
                       error = function(e) { cat("  AN5 fit error:", conditionMessage(e), "\n"); NULL })
            }
            fits <- list(
              `a_primary_diag1e-3` = s_gwas,
              `e_diag1e-3_dropDentist5.45` = if (any(dflag)) an5_fit(R_reg, keep = !dflag) else NULL,
              `f_diag1e-3_L20` = an5_fitL(R_reg, 20L))
            zsets <- list(gwas_z, gwas_z[!dflag], gwas_z)
            S <- list(); C <- list(); P <- list()
            for (k in seq_along(fits)) {
              lab <- names(fits)[k]
              S[[k]] <- an5_summ(fits[[k]], lab, zsets[[k]])
              cc <- an5_coloc(fits[[k]], lab)
              C[[k]] <- cc$s; P[[k]] <- cc$p
            }
            diag_dt <- merge(rbindlist(S, fill = TRUE), rbindlist(C, fill = TRUE), by = "fit", all = TRUE)
            diag_dt <- cbind(gene = gene_symbol, ensembl = gene_id, gwas = gwas_name, chr = chr_num,
                             gwas_n_eff = gwas_n_eff, lambda_s = susie_lambda_s, kr_row, diag_dt)
            fwrite(diag_dt, file.path(an5_dir, paste0(an5_tag, "__fits2.tsv")), sep = "\t")
            pp <- rbindlist(P, fill = TRUE)
            if (nrow(pp)) fwrite(pp, file.path(an5_dir, paste0(an5_tag, "__pairs2.tsv")), sep = "\t")
            # z and LD-implied conditional mean for the region (for plots / later checks)
            cat("  AN5 diagnostics written for", an5_tag, "\n")
          }, error = function(e) cat("  AN5 diagnostic block failed:", conditionMessage(e), "\n"))
          ## ==== end AN5 block ====

