#!/bin/bash
# STAR --readFilesCommand for the S1 mate-order probe: the first 1,000,000 records of a gzipped FASTQ.
zcat "$1" | head -n 4000000 || true
