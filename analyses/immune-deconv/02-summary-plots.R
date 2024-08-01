# Author: Komal S. Rathi
# summarise deconvolution results and create plots

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(rlist)
})

# parse parameters
option_list <- list(
  make_option(c("--deconv_output"), type = "character", 
              help = "deconvolution output from 01-immune.deconv.R (.rds)"),
  make_option(c("--output_dir"), type = "character", 
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
deconv_output <- opt$deconv_output
output_dir <- opt$output_dir

# source plotting theme and heatmap functions
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "immune-deconv", "util", "heatmap_by_group.R"))

# deconvolution method
deconv_output <- readRDS(deconv_output)

# output file
dir.create(output_dir, showWarnings = F, recursive = T)
method <- deconv_output %>%
  pull(method) %>%
  unique()

# annotation colors
annot_colors <- file.path(root_dir, "analyses", "immune-deconv", "util", "colors.yaml")
annot_colors <- list.load(annot_colors)
annot_colors <- lapply(annot_colors, function(x) unlist(x))

# only keep relevant annotation colors  
annot_colors$cohort <- annot_colors$cohort[names(annot_colors$cohort) %in% deconv_output$cohort]
annot_colors$sample_type <- annot_colors$sample_type[names(annot_colors$sample_type) %in% deconv_output$sample_type]

# create heatmap of average immune scores per cell type per cancer and gtex group
output_file <- file.path(output_dir, paste0(method, "_heatmap_by_group.pdf"))
heatmap_by_group(deconv_output = deconv_output, annot_colors = annot_colors,
                 output_file = output_file)
