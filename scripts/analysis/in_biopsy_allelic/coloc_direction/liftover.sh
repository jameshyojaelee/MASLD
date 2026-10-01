#!/bin/bash
# UCSC liftOver in its own module environment: it links OpenSSL 1.1 and MariaDB, while the
# bcftools module used by the caller links OpenSSL 3, so the two cannot share one environment.
set -euo pipefail
source /etc/profile.d/lmod.sh >/dev/null 2>&1 || true
module purge >/dev/null 2>&1
module load Kent_tools/461-GCC-12.2.0 >/dev/null 2>&1
exec liftOver "$@"
