if(!require(optparse)){install.packages('optparse')}
library("optparse")

#process inputs
option_list <- list(
  make_option(
    opt_str = "--input-tsv",
    type = "character",
    dest = "input_tsv",
    help = "tsv RNA expression matrix",
  ),
  make_option(
    opt_str = "--gene-col",
    type = "character",
    dest = "gene_col",
    help = "Name of column with gene symbol"
  ),
  make_option(
    opt_str = "--output-filename",
    type = "character",
    dest = "out",
    help = "Output file name, should fit pattern gene-[expression/counts]-rsem-[fpkm/tpm/expected_count]-collapsed.rds"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
message("Reading in tsv...")
tsv_df = read.table(file = opts$input_tsv, sep = '\t', header = TRUE, check.names = FALSE, row.names=opts$gene_col)
message("Saving as rds...")
saveRDS(tsv_df, file=opts$out)
message("Done!")
