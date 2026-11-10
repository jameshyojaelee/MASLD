# Check the actual canonical/helper functions on a planted stage effect.
source('Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R')
getfun <- function(file, name) {
  for (e in parse(file = file)) {
    if (is.call(e) && identical(e[[1]], as.name('<-')) && identical(e[[2]], as.name(name))) {
      return(eval(e, .GlobalEnv))
    }
  }
  stop('Function not found: ', name)
}
set.seed(20261107)
info <- data.table(batch = rep(1:3, each = 20), inferred_sex = rep(c('F', 'M'), 30),
                   bin = factor(rep(rep(BINS, each = 2), 10), levels = BINS))
K <- paste0('g', 1:40)
RR <- matrix(rnorm(2400), 40, 60) + outer(seq(.2, 2, length.out = 40), as.integer(info$bin))
rownames(RR) <- K; AE <- seq(4, 8, length.out = 40)
ratio_fit <- getfun('Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/code/refit_eisa_adjusted.R', 'ratio_fit')
eisa_fit <- getfun('Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_eisa.R', 'eisa_fit')
obs <- ratio_fit(info$bin)
b <- as.character(info$bin); strat <- paste(info$batch, info$inferred_sex)
for (s in unique(strat)) { j <- which(strat == s); b[j] <- b[j][sample.int(length(j))] }
p <- ratio_fit(b); can_obs <- eisa_fit(info$bin, extra = TRUE); can_perm <- eisa_fit(b, extra = TRUE)
stopifnot(identical(obs, ratio_fit(as.character(info$bin))), max(abs(p$dEmI - obs$dEmI)) > 1e-8,
          max(abs(obs$dEmI - can_obs$dEmI)) < 1e-12, max(abs(p$dEmI - can_perm$dEmI)) < 1e-12,
          all(table(strat, factor(b, levels = BINS)) == table(strat, factor(info$bin, levels = BINS))))
cat('PASS: identity reproduction, actual label perturbation, stratum margins, canonical/helper observed and permuted agreement\n')
