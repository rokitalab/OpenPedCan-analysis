# Author: Komal S. Rathi updated, 2020-07 Kelsey Keith
# Modified By: Arun Boddapati, 2025-01
# script to analyze existing and pre-processed results from xCell and quantiseq.

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggthemes)
  library(ggplot2)
  library(ggpubr) ##Please install this library into docker
  library(ggforce) ##Please install this library into docker

})

# parse parameters
option_list = list(
  make_option(c("--deconv_result_file"), type = "character",
              help = "deconv_method_output(.rds)"),
  make_option(c("--output_dir"), type = "character", 
              help = "output directory"),
  make_option(c ("--filter_cancer_group"),type="character", default = "Medulloblastoma",
              help = "cancer group (e.g., medulloblastoma)"),
  make_option(c("--filter_for_cancer_subtype"),type="character", default = "MB, Group3;MB, Group4;MB, SHH;MB, WNT",
              help = "subtype (e.g., MB,WNT; MB,SHH;MB, Group3)")
              ##           
)

opt = parse_args(OptionParser(option_list = option_list))
deconv_result_file = opt$deconv_result_file
output_dir = opt$output_dir
filter_this_cancer_group = opt$filter_cancer_group #AKB
filter_for_these_cancer_subtype = unlist(strsplit(opt$filter_for_cancer_subtype, ";")) #AKB
print(filter_for_these_cancer_subtype)

deconv_method = gsub( "_output.rds" , "" , basename(opt$deconv_result) )

# output file
dir.create(output_dir, showWarnings = F, recursive = T)
deconv_summary_file = file.path(output_dir, paste0(deconv_method,"_deconv_summary_table.tsv"))

# read deconvulution result file 
deconv_result = readRDS(deconv_result_file)
deconv_summary = deconv_result %>% group_by(cancer_group, cohort, cell_type, molecular_subtype) %>% tally()
write_tsv(deconv_summary, file = deconv_summary_file)

full_output_filtered = deconv_result %>% 
                          filter( cancer_group == filter_this_cancer_group & 
                           molecular_subtype %in% filter_for_these_cancer_subtype & 
                             ! cell_type %in% c("stroma score","immune score","microenvironment score" ) )

full_output_filtered_summary = full_output_filtered %>% group_by(cancer_group, cohort, cell_type, molecular_subtype) %>% tally() %>% data.frame()
filtered_output_summary_file = file.path(output_dir, paste0(deconv_method,"_filtered_output_summary_table.tsv"))
write_tsv(full_output_filtered_summary, file = filtered_output_summary_file)

AllSubtypePairs = combn( unique(full_output_filtered$molecular_subtype) , 2 ,simplify = FALSE)
cohort_groups = unique(full_output_filtered$cohort)

pdf(file.path(output_dir, paste0(deconv_method,"_immuneCell_distributions_byMolecularSubtype.pdf")), width = 20, height = 25)

for (cohort_name in cohort_groups) {
  # Filter data for the current cohort
  # Generate plot
      
  p =  full_output_filtered %>% 
                              filter( cohort == cohort_name ) %>%
                                    ggplot( aes(x = molecular_subtype, y = fraction, fill = molecular_subtype)) + 
                                    geom_boxplot(notch = TRUE) + 
                                    facet_wrap(~ cell_type) +
                                    theme_bw() +
                                    theme(
                                      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                      legend.position = "top"
                                    ) +
                                    stat_compare_means(method = "wilcox", comparisons = AllSubtypePairs) +  
                                    labs(fill = "Molecular Subtype", title = paste0(filter_this_cancer_group , "\n" , cohort_name))
  
  # Print plot to PDF
  print(p)
}

# Close the PDF device
dev.off()


