library(tidyverse)
library(dplyr)

<<<<<<< HEAD
=======
#v27_eh_df <- read_tsv('results/ensembl_gene_symbol_gtf_gencode_v27.tsv')
#v28_eh_df <- read_tsv('results/ensembl_gene_symbol_gtf_gencode_v28.tsv')
#v36_eh_df <- read_tsv('results/ensembl_gene_symbol_gtf_gencode_v36.tsv')
>>>>>>> 142399b1dd30453568dd8509917cc11477fc7be7
v39_eh_df <- read_tsv('results/ensembl_gene_symbol_gtf_gencode_v39.tsv')

open_ped_can_eh_df <- read_tsv('input/OpenPedCan-v12-ensg-hugo-pmtl-mapping.tsv') %>%
  dplyr::select(gene_symbol, ensg_id) %>%
  dplyr::rename(ensembl = ensg_id)

<<<<<<< HEAD
=======
#merged_eh_df <- distinct(bind_rows(v27_eh_df, v28_eh_df, v36_eh_df, v39_eh_df, open_ped_can_v7_eh_df))
>>>>>>> 142399b1dd30453568dd8509917cc11477fc7be7
merged_eh_df <- distinct(bind_rows(v39_eh_df, open_ped_can_eh_df))

write_tsv(merged_eh_df, 'results/gencode_ensg_symbol_map_merged.tsv')
