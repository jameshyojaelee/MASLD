#!/usr/bin/env Rscript
# Arm A2: the published ordering WITHOUT ComBat.
#
# Identical to A1 in every other respect, so the A1/A2 difference isolates one
# thing: what batch correction does to the axis.
#
# The reason this arm exists. Kamzolas applied ComBat across cohorts that each
# spanned the disease spectrum, so removing cohort could not remove much disease.
# Our substrate is not like that: GSE126848 (53) is entirely unstaged control and
# obese donors, and GSE213621 (361) is entirely staged disease. Here cohort partly
# IS disease, so ComBat with sex but not disease in the model can remove real
# signal. A2 keeps cohort out of the transform and handles it as a model covariate
# downstream instead.
#
# This is a wrapper rather than a copy so the two arms can never drift apart.

Sys.setenv(CAB_ARM = "A2", CAB_MATRIX = "E_qn")
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/02a_arm_kamzolas_exact.R"))
