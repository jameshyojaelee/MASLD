#!/bin/bash
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/downloads_untrusted
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
curl -sS --fail -o $D/height_GCST90245992/GCST90245992_buildGRCh37.tsv $B/GCST90245001-GCST90246000/GCST90245992/GCST90245992_buildGRCh37.tsv
curl -sS --fail -o $D/height_GCST90245992/md5sum.txt $B/GCST90245001-GCST90246000/GCST90245992/md5sum.txt
curl -sS --fail -o $D/height_GCST90245992/meta.yaml $B/GCST90245001-GCST90246000/GCST90245992/GCST90245992_buildGRCh37.tsv-meta.yaml
curl -sS --fail -o $D/fi_GCST90002238/GCST90002238_buildGRCh37.tsv.gz $B/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz
curl -sS --fail -o $D/fi_GCST90002238/md5sum.txt $B/GCST90002001-GCST90003000/GCST90002238/md5sum.txt
curl -sS --fail -o $D/fi_GCST90002238/README.txt $B/GCST90002001-GCST90003000/GCST90002238/README.txt
curl -sS --fail -o $D/dpw_GCST007461/DrinksPerWeek.txt.gz $B/GCST007001-GCST008000/GCST007461/DrinksPerWeek.txt.gz
curl -sS --fail -o $D/dpw_GCST007461/README $B/GCST007001-GCST008000/GCST007461/README
curl -sS --fail -o $D/tg_GCST90239664/GCST90239664_buildGRCh37.tsv $B/GCST90239001-GCST90240000/GCST90239664/GCST90239664_buildGRCh37.tsv
curl -sS --fail -o $D/tg_GCST90239664/md5sum.txt $B/GCST90239001-GCST90240000/GCST90239664/md5sum.txt
curl -sS --fail -o $D/tg_GCST90239664/README.txt $B/GCST90239001-GCST90240000/GCST90239664/README.txt
for d in height_GCST90245992 fi_GCST90002238 tg_GCST90239664; do (cd $D/$d && grep -v -e meta.yaml -e README -e "md5sum.txt" md5sum.txt | md5sum -c -); done
md5sum $D/dpw_GCST007461/DrinksPerWeek.txt.gz
du -sh $D/*
echo DOWNLOADS DONE
