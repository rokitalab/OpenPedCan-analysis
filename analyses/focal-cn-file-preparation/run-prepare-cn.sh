#!/bin/bash
# C. Bethell and C. Savonen for CCDL 2019, J. Rokita for D3b 2023
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# Run original files - will not by default
RUN_ORIGINAL=${RUN_ORIGINAL:-0}

# Run testing files for circle CI - will not by default
IS_CI=${OPENPBTA_TESTING:-0}

RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
results_dir=../../analyses/focal-cn-file-preparation/results
gtf_file=${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz

if [[ "$RUN_FOR_SUBTYPING" -eq "1" ]]
then
  histologies_file=${data_dir}/histologies-base.tsv
else
  histologies_file=${data_dir}/histologies.tsv
fi

echo "add ploidy to consensus"
# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

echo "running bedtools coverage for each sample bed"
# Run snakemake script implementing `bedtools coverage` for each sample bed file in
# `scratch/cytoband_status` -- these files are generated in
# `02-add-ploidy-consensus.Rmd`
# currently runs 10 jobs in parallel, which should be fine for most implementations
snakemake -j 10 --snakefile run-bedtools.snakemake

echo "Determine the dominant status for each chromosome arm and compare with GISTIC's arm status calls"
Rscript --vanilla -e "rmarkdown::render('03-add-cytoband-status-consensus.Rmd', clean = TRUE)"

echo "Run annotation step for consensus file"
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${results_dir}/consensus_seg_with_status.tsv \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "consensus_seg_annotated_cn" \
--seg

echo "Define most focal units of recurrent CNVs"
Rscript --vanilla -e "rmarkdown::render('05-define-most-focal-cn-units.Rmd', clean = TRUE)"

echo "Define the recurrent calls"
Rscript --vanilla -e "rmarkdown::render('06-find-recurrent-calls.Rmd', clean = TRUE)"


# if we want to process the CNV data from the original callers
# (e.g., CNVkit, ControlFreeC)
if [ "$RUN_ORIGINAL" -gt "0" ]
then
  
  echo "prep CNVkit WXS data"
  # Prep the CNVkit data
  Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"
  
  echo "annotate CNVkit WXS CNVs"
  # Run annotation step for CNVkit WXS only
  Rscript --vanilla 04-prepare-cn-file.R \
  --cnv_file ${results_dir}/cnvkit_with_status.tsv \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "cnvkit_annotated_cn" \
  --seg \
  --runWXSonly
  
  echo "annotate FREEC WXS CNVs"
  # Run annotation step for ControlFreeC
  Rscript --vanilla 04-prepare-cn-file.R \
  --cnv_file ${data_dir}/cnv-controlfreec.tsv.gz \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "controlfreec_annotated_cn" \
  --controlfreec \
  --runWXSonly
  
  echo "annotate CNVkit tumor only CNVs"
  # Run annotation step for ControlFreeC tumor only
  Rscript --vanilla 04-prepare-cn-file.R \
  --cnv_file ${data_dir}/cnv-controlfreec-tumor-only.tsv.gz \
  --gtf_file $gtf_file \
  --metadata $histologies_file \
  --filename_lead "freec-tumor-only_annotated_cn" \
  --controlfreec

results_dir=./results

  echo "merge annotated files with cnvkit WXS"
  # Run merging for all annotated files 
  Rscript --vanilla 07-consensus-annotated-merge.R \
  --cnvkit_auto ${results_dir}/cnvkit_annotated_cn_wxs_autosomes.tsv.gz \
  --cnvkit_x_and_y ${results_dir}/cnvkit_annotated_cn_wxs_x_and_y.tsv.gz \
  --consensus_auto ${results_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz \
  --consensus_x_and_y ${results_dir}/consensus_seg_annotated_cn_x_and_y.tsv.gz \
  --cnv_tumor_auto ${results_dir}/freec-tumor-only_annotated_cn_autosomes.tsv.gz \
  --cnv_tumor_x_and_y ${results_dir}/freec-tumor-only_annotated_cn_x_and_y.tsv.gz \
  --cnvkitWXS TRUE \
  --outdir ${results_dir}
  
  echo "merge annotated files with freec WXS"
  # Run merging for all annotated files 
  Rscript --vanilla 07-consensus-annotated-merge.R \
  --freec_auto ${results_dir}/controlfreec_annotated_cn_wxs_autosomes.tsv.gz \
  --freec_x_and_y ${results_dir}/controlfreec_annotated_cn_wxs_x_and_y.tsv.gz \
  --consensus_auto ${results_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz \
  --consensus_x_and_y ${results_dir}/consensus_seg_annotated_cn_x_and_y.tsv.gz \
  --cnv_tumor_auto ${results_dir}/freec-tumor-only_annotated_cn_autosomes.tsv.gz \
  --cnv_tumor_x_and_y ${results_dir}/freec-tumor-only_annotated_cn_x_and_y.tsv.gz \
  --cnvkitWXS FALSE \
  --outdir ${results_dir}

fi
