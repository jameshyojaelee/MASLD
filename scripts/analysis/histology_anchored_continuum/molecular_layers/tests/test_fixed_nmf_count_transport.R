suppressPackageStartupMessages(library(edgeR))
source("scripts/analysis/histology_anchored_continuum/molecular_layers/fixed_nmf_count_transport.R")
annotation <- data.frame(gene_id = c("G1", "A1", "G2", "G3", "G4", "G5"),
  counting_gene_id = c("G1", "G1", "G2", "G3", "", "G5"),
  disposition = c("primary_gene", "equivalent", "primary_gene", "unresolved", "excluded_reference", "duplicated"))
mask <- fixed_nmf_counting_mask(c("G1", "A1", "G2", "G3", "G4", "G5", "absent"), annotation)
stopifnot(identical(mask$retained, c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)))
# An ALT outside the original basis must also mask its primary counting target.
stopifnot(!fixed_nmf_counting_mask("G1", annotation)$retained)
set.seed(20260923)
reference <- setNames(rpois(200, 100), paste0("G", 1:200))
counts <- cbind(S1 = reference * 2, S2 = reference * 3)
sizes <- setNames(colSums(counts), colnames(counts))
nf <- fixed_reference_tmm(counts, sizes, reference, sum(reference), .9)
single <- fixed_reference_tmm(counts[, "S1", drop = FALSE], sizes["S1"], reference, sum(reference), .9)
stopifnot(abs(nf[["S1"]] - single[["S1"]]) < 1e-12, all(abs(nf - .9) < 1e-12))
state <- list(row_minimum = setNames(rep(0, 200), names(reference)),
  cohort_batch_subtraction = matrix(0, 200, 1, dimnames = list(names(reference), "C1")))
result <- fixed_nmf_transport_expression(counts, sizes, nf, state, c("C1", "C1"))
stopifnot(all(is.finite(result)), all(result >= 0), identical(dimnames(result), dimnames(counts)))
unknown <- tryCatch(fixed_nmf_transport_expression(counts, sizes, nf, state, c("unknown", "C1")), error = identity)
stopifnot(inherits(unknown, "error"))
cat("PASS: counting-object mask, ALT contributors outside basis, fixed-reference sample independence, finite transform and unknown-cohort rejection\n")
# Reconstruction criteria sit exactly on their boundaries and are not relaxed.
stopifnot(fixed_nmf_axis_eligible(.99, .9, .8, 4L, 4L), !fixed_nmf_axis_eligible(.9899, .95, .9, 4L, 4L),
  !fixed_nmf_axis_eligible(.995, .8999, .9, 6L, 6L), !fixed_nmf_axis_eligible(.995, .95, .7999, 6L, 6L),
  !fixed_nmf_axis_eligible(.995, .95, .9, 5L, 6L))
checks <- data.frame(k = c(4L, 4L, 6L, 6L, 4L), program = c("P1", "P2", "P5", "P6", "P2"),
  variant = c("full_stage", "full_stage", "full_stage", "full_stage", "signature_excluded_continuum"),
  transport_reconstruction_eligible = c(TRUE, FALSE, TRUE, FALSE, TRUE))
stopifnot(identical(fixed_nmf_unavailable_axes(checks, "full_stage"), c("k4-P2", "k6-P6")),
  identical(fixed_nmf_unavailable_axes(checks, "signature_excluded_continuum"), character()))
# NNLS recovers a known nonnegative mixture and rejects mismatched feature rows.
set.seed(20260923)
w <- matrix(rexp(60 * 3), 60, 3, dimnames = list(paste0("G", 1:60), NULL))
h_true <- matrix(c(1, 0, 2, 0, 3, 1, 0, 0, 4, 1, 1, 0), 3, 4)
x <- w %*% h_true; dimnames(x) <- list(rownames(w), paste0("S", 1:4))
h <- fixed_nmf_nnls_scores(w, x)
stopifnot(identical(rownames(h), c("P1", "P2", "P3")), identical(colnames(h), colnames(x)), max(abs(h - h_true)) < 1e-8)
bad <- tryCatch(fixed_nmf_nnls_scores(w, x[rev(rownames(x)), ]), error = identity)
stopifnot(inherits(bad, "error"))
# Unavailable axes become NA for every sample; every other row is untouched.
masked <- fixed_nmf_mask_unavailable(h, 3L, c("k3-P2", "k6-P6"))
stopifnot(all(is.na(masked["P2", ])), identical(masked[c("P1", "P3"), ], h[c("P1", "P3"), ]),
  identical(fixed_nmf_mask_unavailable(h, 3L, character()), h))
wrong_k <- tryCatch(fixed_nmf_mask_unavailable(h, 4L, "k4-P2"), error = identity)
stopifnot(inherits(wrong_k, "error"))
cat("PASS: axis criteria boundaries, unavailable-axis list, NNLS recovery, row mismatch rejection, whole-axis NA masking\n")
