#!/usr/bin/env Rscript
# Arm A2b: the published ordering, with disease PROTECTED in the ComBat model.
#
# A1 used ComBat(mod = ~ sex), as published. In this substrate cohort is partly
# disease, and that ComBat drove PC1's cohort eta-squared to zero, so it removed
# every between-cohort difference including any that was real. A2b keeps
# fibrosis stage in `mod` so ComBat cannot remove it. The A1/A2b difference
# therefore isolates one thing: whether the Q1 null was an over-correction
# artifact.
#
# A wrapper rather than a copy, exactly like A2, so the three arms can never
# drift apart. Requires 01b to have written arms/E_cbd.rds first.
#
# READ 01b's HEADER BEFORE READING ANY A2b NUMBER. This arm is post-hoc, cannot
# pass the sealed Q1 gate, and is circular on its own: protecting stage
# preserves stage-associated variance by construction. The interpretable
# quantity is A2b minus the permuted-protection null that 02i builds, not A2b
# alone.

Sys.setenv(CAB_ARM = "A2b", CAB_MATRIX = "E_cbd")
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/02a_arm_kamzolas_exact.R"))
