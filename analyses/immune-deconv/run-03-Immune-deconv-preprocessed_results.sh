#!/bin/bash
# Module author: Komal S. Rathi, updated Kelsey Keith
# Modified by: Arun Boddapati
# Modified on: 2025-04

# This script plots the fraction of immune cell types across cancer molecular subtypes
# for the pre-processed immune deconv results (xCell and quantiSeq)

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
echo "Plotting fractions for cancer molecular subtype for quantiseq deconv method"
Rscript --vanilla 03-Immune-deconv-preprocessed_results.R \
--deconv_result 'results/quantiseq_output.rds' \
--output_dir 'results_custom_deconv_subtyping'

### quanTIseq
#generate deconvolution output
echo "Plotting fractions for cancer molecular subtype for xcell deconv method"
Rscript --vanilla 03-Immune-deconv-preprocessed_results.R \
--deconv_result 'results/xcell_output.rds' \
--output_dir 'results_custom_deconv_subtyping'