#!/bin/bash
# Module author: Komal S. Rathi, updated Kelsey Keith
# 2022-07
# Modified from original author to add on immune enrichment analysis script (immune-enrichment.R) - Audrey Crowther 2025-01

# This script runs the steps for immune deconvolution using xCell and quanTIseq. 

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

### xCell
# generate deconvolution output
echo "Deconvolution xCell..."
Rscript --vanilla 01-immune-deconv.R \
--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
--clin_file '../../data/histologies.tsv' \
--deconv_method 'xcell' \
--output_dir 'results'

### quanTIseq
# generate deconvolution output
#echo "Deconvolution quanTIseq..."
#Rscript --vanilla 01-immune-deconv.R \
#--expr_mat '../../data/gene-expression-rsem-tpm-collapsed.rds' \
#--clin_file '../../data/histologies.tsv' \
#--deconv_method 'quantiseq' \
#--output_dir 'results'

### Medullablastoma immune enrichment from xCell results
echo "Medullablastoma Subtypes - Immune Enrichments"
Rscript --vanilla immune-enrichment.R \