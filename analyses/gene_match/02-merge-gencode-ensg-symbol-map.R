library(tidyverse)
library(dplyr)

v39_eh_df <- read_tsv('results/ensembl_gene_symbol_gtf_gencode_v39.tsv')

open_ped_can_eh_df <- read_tsv('input/OpenPedCan-v12-ensg-hugo-pmtl-mapping.tsv') %>%
  dplyr::select(gene_symbol, ensg_id) %>%
  dplyr::rename(ensembl = ensg_id)

merged_eh_df <- distinct(bind_rows(v39_eh_df, open_ped_can_eh_df))

write_tsv(merged_eh_df, 'results/gencode_ensg_symbol_map_merged.tsv')
