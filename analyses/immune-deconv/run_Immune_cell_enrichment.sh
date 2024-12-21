#!/bin/bash

# Define input arguments
Xcell_output="/path/to/xcell_output.rds"  # Replace with the actual path
clin_file="/path/to/histologies.tsv"      # Replace with the actual path
output_dir="/path/to/output_directory"    # Replace with the actual path

# Run the R script
Rscript 01-Immune_cell_enrichment_in_MB_subtypes.R \
  --Xcell_output $Xcell_output \
  --clin_file $clin_file \
  --output_dir $output_dir
