# Author: Prawat ###
# Immune Cell Enrichment Analysis in Medulloblastoma Cancer Subtypes ####
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
xcell_output_path <- args$Xcell_output
clin_file_path <- args$clin_file
output_dir <- args$output_dir

# Load xcell deconvolution results and clinical data
xcell_output <- readRDS(xcell_output_path)
clin_file <- read_tsv(clin_file_path, guess_max = 10000)

# Filter Data for Medulloblastoma Samples
medulloblastoma_samples <- clin_file %>%
  filter(pathology_diagnosis == "Medulloblastoma" & experimental_strategy == "RNA-Seq")

filtered_xcell_output_MB <- xcell_output %>%
  filter(Kids_First_Biospecimen_ID %in% medulloblastoma_samples$Kids_First_Biospecimen_ID)

# Ensure 'fraction' column is numeric
filtered_xcell_output_MB$fraction <- as.numeric(filtered_xcell_output_MB$fraction)

# Identify unique cell types and molecular subtypes
unique_cell_types <- filtered_xcell_output_MB %>%
  distinct(cell_type) %>%
  pull(cell_type)

unique_molecular_subtypes <- filtered_xcell_output_MB %>%
  distinct(molecular_subtype) %>%
  pull(molecular_subtype)

cat("Unique Cell Types:", paste(unique_cell_types, collapse = ", "), "\n")
cat("Unique Molecular Subtypes:", paste(unique_molecular_subtypes, collapse = ", "), "\n")

# Add cell grouping to simplify 
filtered_xcell_output_MB <- filtered_xcell_output_MB %>%
  mutate(grouped_cell_type = case_when(
    str_detect(cell_type, "T cell|NK") ~ "Lymphoid cells",
    str_detect(cell_type, "B cell") ~ "Lymphoid cells",
    str_detect(cell_type, "Monocyte|Macrophage|Neutrophil|Mast cell") ~ "Myeloid cells",
    str_detect(cell_type, "Progenitor|Stem cell") ~ "Progenitors",
    str_detect(cell_type, "Fibroblast") ~ "Fibroblasts",
    str_detect(cell_type, "Stroma|Microenvironment") ~ "Microenvironment",
    TRUE ~ "Other"
  ))

# Perform unpaired t-tests for each cell type across subtypes
t_test_results <- filtered_xcell_output_MB %>%
  group_by(cell_type) %>%
  summarize(
    t_test = list(combn(unique(molecular_subtype), 2, function(pair) {
      group1 <- filter(filtered_xcell_output_MB, cell_type == first(cell_type), molecular_subtype == pair[1])
      group2 <- filter(filtered_xcell_output_MB, cell_type == first(cell_type), molecular_subtype == pair[2])
      t_test <- t.test(group1$fraction, group2$fraction, var.equal = TRUE)
      data.frame(
        Comparison = paste(pair[1], "vs", pair[2]),
        P_Value = t_test$p.value,
        Mean_Group1 = mean(group1$fraction, na.rm = TRUE),
        Mean_Group2 = mean(group2$fraction, na.rm = TRUE)
      )
    }, simplify = FALSE)),
    .groups = "drop"
  ) %>%
  unnest(cols = c(t_test))

# Adjust p-Values using Benjamini-Hochberg Method
t_test_results <- t_test_results %>%
  mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH"))

# Save t-test results 
write_csv(t_test_results, file.path(output_dir, "unpaired_t_test_results_MB.csv"))

# Generate boxplots with significance comparisons for each cell type
if (nrow(t_test_results) > 0) {
  t_test_results <- t_test_results %>%
    mutate(Significance = case_when(
      Adjusted_P_Value < 0.001 ~ "***",
      Adjusted_P_Value < 0.01 ~ "**",
      Adjusted_P_Value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  # Create boxplots with Immune cell distribution among molecular subtypes
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
  
  # Save the Plot
  ggsave(file.path(output_dir, "Enriched Immune cell types in Medulloblastoma subtypes.pdf"), plot, width = 24, height = 20)
}

# Completion Message
cat("Analysis complete. Results and plots saved to:", output_dir, "\n")
