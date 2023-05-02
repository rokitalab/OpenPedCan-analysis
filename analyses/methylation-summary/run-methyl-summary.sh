#!/bin/bash
# PediatricOpenTargets 2022
# Eric Wafula
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up paths to input and output directories for the analysis
absolute_path=$(cd "$(dirname "$0")"; pwd -P)
data_dir="${absolute_path}/../../data"
#input_dir="input"
results_dir="results"

printf '\nStarting Analysis...\n'

# Create tmp/ directory for R scripts
mkdir -p -m777 ./tmp

###################### Calculate gene-level tmp median expression ##########################
TMP=./tmp TMPDIR=./tmp Rscript --vanilla 01-calculate-tpm-medians.R  \
--histologies ${data_dir}/histologies.tsv \
--rnaseq_matrix ${data_dir}/gene-expression-rsem-tpm-collapsed.rds \
--methyl_probe_annot ${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz \
--methyl_independent_samples ${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
--rnaseq_independent_samples ${data_dir}/independent-specimens.rnaseqpanel.primary.eachcohort.tsv \
--exp_values gene

###################### Calculate isoform-level tmp median expression ##########################
TMP=./tmp TMPDIR=./tmp Rscript --vanilla 01-calculate-tpm-medians.R  \
--histologies ${data_dir}/histologies.tsv \
--rnaseq_matrix ${data_dir}/rna-isoform-expression-rsem-tpm.rds \
--methyl_probe_annot ${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz \
--methyl_independent_samples ${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
--rnaseq_independent_samples ${data_dir}/independent-specimens.rnaseqpanel.primary.eachcohort.tsv \
--exp_values isoform

###################### Calculate probe-level beta quantiles #################################
TMP=./tmp TMPDIR=./tmp Rscript --vanilla 02-calculate-methly-quantiles.R \
--histologies ${data_dir}/histologies.tsv \
--methyl_matrix ${data_dir}/methyl-beta-values.rds \
--methyl_probe_annot ${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz \
--independent_samples ${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
--methyl_values beta

###################### Calculate correlations between beta and gene tpm values ###############
python3 03-methyl-tpm-correlation.py \
--methyl_values beta \
--exp_values gene \
${data_dir}/histologies.tsv \
${data_dir}/independent-specimens.rnaseqpanel.primary.eachcohort.tsv \
${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
${data_dir}/methyl-beta-values.rds \
${data_dir}/gene-expression-rsem-tpm-collapsed.rds \
${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz

###################### Calculate correlations between beta and isoform tpm values ###########
python3 03-methyl-tpm-correlation.py \
--methyl_values beta \
--exp_values isoform \
${data_dir}/histologies.tsv \
${data_dir}/independent-specimens.rnaseqpanel.primary.eachcohort.tsv \
${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
${data_dir}/methyl-beta-values.rds \
${data_dir}/rna-isoform-expression-rsem-tpm.rds \
${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz

############################ Create expression tpm transcript representation ################
python3 04-tpm-transcript-representation.py \
${data_dir}/histologies.tsv \
${data_dir}/independent-specimens.rnaseqpanel.primary.eachcohort.tsv \
${data_dir}/independent-specimens.methyl.primary.eachcohort.tsv \
${data_dir}/gene-expression-rsem-tpm-collapsed.rds \
${data_dir}/rna-isoform-expression-rsem-tpm.rds \
${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz

############################ Create gene-level beta methylation summary table ################
TMP=./tmp TMPDIR=./tmp Rscript --vanilla 05-create-methyl-summary-table.R \
--methyl_tpm_corr ${results_dir}/gene-methyl-probe-beta-tpm-correlations.tsv.gz \
--methyl_probe_qtiles ${results_dir}/methyl-probe-beta-quantiles.tsv.gz \
--methyl_probe_annot ${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz \
--rnaseq_tpm_medians ${results_dir}/gene-median-tpm-expression.tsv.gz \
--efo_mondo_annot ${data_dir}/efo-mondo-map.tsv \
--exp_values gene \
--methyl_values beta

############################ Create isoform-level beta methylation summary table ##############
TMP=./tmp TMPDIR=./tmp Rscript --vanilla 05-create-methyl-summary-table.R \
--methyl_tpm_corr ${results_dir}/isoform-methyl-probe-beta-tpm-correlations.tsv.gz \
--methyl_probe_qtiles ${results_dir}/methyl-probe-beta-quantiles.tsv.gz \
--methyl_probe_annot ${data_dir}/infinium.gencode.v39.probe.annotations.tsv.gz \
--rnaseq_tpm_medians ${results_dir}/isoform-median-tpm-expression.tsv.gz \
--tpm_transcript_rep ${results_dir}/methyl-tpm-transcript-representation.tsv.gz \
--efo_mondo_annot ${data_dir}/efo-mondo-map.tsv \
--exp_values isoform \
--methyl_values beta


########################### Convert methytion summary table from TSV to (JSONL) ###############
python3 06-methly-summary-tsv2jsonl.py \
--methyl_values beta \
${results_dir}/gene-methyl-beta-values-summary.tsv.gz \
${results_dir}/isoform-methyl-beta-values-summary.tsv.gz

# remove tmp/ directory
rm -rf ./tmp

printf '\nAnalysis Done...\n'
