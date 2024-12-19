#!/bin/bash
Rscript Immune_cell_enrichment_MB_subtypes.R \
  --expr_mat /path/to/expression_matrix.rds \
  --clin_file /path/to/histologies.tsv \
  --deconv_method xcell \
  --output_dir /path/to/output_directory
