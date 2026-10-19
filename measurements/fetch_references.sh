#!/bin/bash
# Fetch published reference tables (large files not tracked in git).
set -e
D=$(dirname $(readlink -f $0))
cd $D/ccqelike_3d_me_fhc/reference
for f in README MINERvA_TripleDiffQELike_ptpzsumtp.root \
         MINERvA_TripleDiffQELike_sumtpptpz_data_cross_section.csv \
         MINERvA_TripleDiffQELike_sumtpptpz_covariances_underoverbinsremoved.csv; do
  [ -f $f ] || curl -sfL https://arxiv.org/src/2203.08022v1/anc/$f -o $f
done
