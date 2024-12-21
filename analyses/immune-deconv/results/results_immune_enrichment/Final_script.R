# Author: Prawat ###
# Comprehensive Immune Cell Enrichment Analysis in Medulloblastoma Cancer Subtypes
# Date: 2024-12-20

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggpubr)
  library(optparse)
})

# Parse input arguments
options_list <- list(
  make_option(c("-e", "--Xcell_output"), type = "character",
              help = "Path to the xcell_output (.rds)"),
  make_option(c("-c", "--clin_file"), type = "character",
              help = "Path to the histologies metadata file (.tsv)"),
  make_option(c("-o", "--output_dir"), type = "character",
              help = "Directory to save output files")
)

args <- parse_args(OptionParser(option_list = options_list))

# Assign parameters
#xcell_output_path <- args$Xcell_output
xcell_output_path<-"/home/rstudio/OpenPedCan-Project/analyses/immune-deconv/results/xcell_output.rds"
#clin_file_path <- args$clin_file
clin_file_path<-"/home/rstudio/OpenPedCan-Project/data/histologies.tsv"
#output_dir <- args$output_dir
output_dir<-"/home/rstudio/OpenPedCan-Project/analyses/immune-deconv/results/results_immune_enrichment/Final_results/"
dir.create(output_dir, showWarnings = FALSE)

# Load xcell deconvolution results and clinical data
xcell_output <- readRDS(xcell_output_path)
clin_file <- read_tsv(clin_file_path, guess_max = 10000)

# Filter clinical data for Medulloblastoma Samples
medulloblastoma_samples <- clin_file %>%
  filter(pathology_diagnosis == "Medulloblastoma" & experimental_strategy == "RNA-Seq")

#Filter enrichment data for Medulloblastoma samples#########
filtered_xcell_output_MB <- xcell_output %>%
  filter(Kids_First_Biospecimen_ID %in% medulloblastoma_samples$Kids_First_Biospecimen_ID)

# Ensure 'fraction' column is numeric
filtered_xcell_output_MB$fraction <- as.numeric(filtered_xcell_output_MB$fraction)

# Add cell grouping
filtered_xcell_output_MB <- filtered_xcell_output_MB %>%
  mutate(grouped_cell_type = case_when(
    str_detect(cell_type, "T cell|NK") ~ "Lymphoid cells",
    str_detect(cell_type, "B cell") ~ "Lymphoid cells",
    str_detect(cell_type, "Monocyte|Macrophage|Neutrophil|Mast cell") ~ "Myeloid cells",
    str_detect(cell_type, "fibroblast") ~ "Fibroblasts",
    str_detect(cell_type, "progenitor|stem cell") ~ "Progenitors",
    str_detect(cell_type, "stroma|microenvironment") ~ "Microenvironment",
    TRUE ~ "Other"
  ))

# 1. Generate Boxplots with Significance Comparisons
cat("Generating boxplots with significance comparisons...\n")
# Perform unpaired t-tests for each cell type across subtypes
t_test_results <- filtered_xcell_output_MB %>%
  group_by(cell_type) %>%
  summarize(
    t_test = list(combn(unique(molecular_subtype), 2, function(pair) {
      group1 <- filter(filtered_xcell_output_MB, cell_type == first(cell_type), molecular_subtype == pair[1])
      group2 <- filter(filtered_xcell_output_MB, cell_type == first(cell_type), molecular_subtype == pair[2])
      
      # Perform t-test
      if (nrow(group1) > 1 && nrow(group2) > 1) {
        t_test <- t.test(group1$fraction, group2$fraction, var.equal = TRUE)
        return(data.frame(
          Comparison = paste(pair[1], "vs", pair[2]),
          P_Value = t_test$p.value,
          Mean_Group1 = mean(group1$fraction, na.rm = TRUE),
          Mean_Group2 = mean(group2$fraction, na.rm = TRUE)
        ))
      } else {
        return(data.frame(
          Comparison = paste(pair[1], "vs", pair[2]),
          P_Value = NA,
          Mean_Group1 = NA,
          Mean_Group2 = NA
        ))
      }
    }, simplify = FALSE)),
    .groups = "drop"
  ) 

# Unnest the `t_test` column
t_test_results <- t_test_results %>%
  unnest(cols = c(t_test))

# Filter out rows with NA P_Value
t_test_results <- t_test_results %>%
  filter(!is.na(P_Value))

# Adjust P-Values using Benjamini-Hochberg Method
t_test_results <- t_test_results %>%
  mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH"))


write_csv(t_test_results, file.path(output_dir, "unpaired_t_test_results_MB.csv"))
#############################################################################################################
######### Plot the box plots for visualizing the enrichment of immune cells ####
plot <- ggplot(filtered_xcell_output_MB, aes(x = molecular_subtype, y = fraction, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(), alpha = 0.5, size = 1) +
  stat_compare_means(
    comparisons = combn(unique(filtered_xcell_output_MB$molecular_subtype), 2, simplify = FALSE),
    label = "p.signif",
    method = "t.test",
    paired = FALSE
  ) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Immune Cell Enrichment by Molecular Subtype in Medulloblastoma",
    x = "Molecular Subtype",
    y = "XCell Enrichment Score",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "immune_cell_enrichment_boxplot.pdf"), plot, width = 24, height = 20)

# 2. Generate Heatmap for Average xCell Scores Across Molecular Subtypes
cat("Generating heatmap for average xCell scores...\n")
aggregated_data <- filtered_xcell_output_MB %>%
  group_by(molecular_subtype, cell_type) %>%
  summarize(avg_fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

heatmap_data <- aggregated_data %>%
  pivot_wider(
    names_from = cell_type,
    values_from = avg_fraction,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "molecular_subtype")

cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}
heatmap_data_norm <- t(apply(heatmap_data, 1, cal_z_score))

pdf(file.path(output_dir, "average_xcell_scores_heatmap.pdf"), width = 14, height = 10)
pheatmap(
  heatmap_data_norm,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Average xCell Scores Across Molecular Subtypes"
)
dev.off()

# Completion Message
cat("Analysis complete. Results saved to:", output_dir, "\n")
