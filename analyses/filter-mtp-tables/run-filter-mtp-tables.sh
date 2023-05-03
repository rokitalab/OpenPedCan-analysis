#!/bin/bash
# PediatricOpenTargets 2021
# Eric
set -e
set -o pipefail

printf "Start filtering mtp tables...\n\n"

# Set mtp tables commit data location on s3 bucket to download files 
openpedcan_url="https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v12/mtp-tables/commit"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


# Set directory in scratch to download mtp commit files from s3 bucket
scratch_dir="../../scratch"
mtp_commit_dir="$scratch_dir/mtp-commit"
mkdir -p $mtp_commit_dir
mtp_filtered_dir="$scratch_dir/mtp-filtered"
mkdir -p $mtp_filtered_dir

# Download mtp tsv tables to be filtered
printf "\nDownloading TSV files to filter..."
cd $mtp_commit_dir

wget $openpedcan_url/md5sum.txt
grep '.tsv.gz' < md5sum.txt > md5sum_filter.txt
FILES=(`tr -s ' ' < md5sum_filter.txt | cut -d ' ' -f 2`)
for file in "${FILES[@]}"
do
  wget $openpedcan_url/$file
done

# Check the md5s for downloaded files
printf "\nChecking MD5 hashes..."
md5sum -c md5sum_filter.txt

# Set up paths results directory
cd $script_directory
results_path="results"

###################### Filter MTP tables for current Gencode ###################
printf '\nFiltering MTP  tables for current OT and targets and diseases...'
cd $script_directory

# Create tmp/ directory for R scripts
mkdir -p -m777 ./tmp

TMP=./tmp TMPDIR=./tmp Rscript -e "rmarkdown::render('01-filter-mtp-tables-for-current-gencode.Rmd', \
  clean = TRUE)"

# remove tmp/ directory
rm -rf ./tmp

###################### Convert TSV to JSON Lines (JSONL) ######################
printf '\nConvert TSV files to JSONL files...\n'
cd $mtp_filtered_dir

FILES=(`ls *.tsv.gz`)
for file in "${FILES[@]}"
do
  python3 $script_directory/02-mtp-tsv2jsonl.py $file
done

printf "Done filtering mtp tables...\n\n"