# # Author: prawat ###
# Identifying immune cell enrichment in Medulloblastoma cancer subtypes ####

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(immunedeconv)
  library(optparse)
  library(ggalluvial)
})

require(optparse)  # Ensure optparse is loaded

# Parse input arguments
options_list <- list(
  make_option(c("-e", "--expr_mat"), type = "character",
              help = "Path to the gene expression matrix (.rds)"),
  make_option(c("-c", "--clin_file"), type = "character",
              help = "Path to the histologies metadata file (.tsv)"),
  make_option(c("-m", "--deconv_method"), type = "character",
              help = "Method for immune deconvolution (xcell or quantiseq)"),
  make_option(c("-o", "--output_dir"), type = "character", 
              help = "Directory to save output files")
)

# Parse the arguments
args <- parse_args(OptionParser(option_list = options_list))

# Assign parameters to variables
expr_mat_path <- args
clin_file_path <- args
deconv_method <- tolower(args)
output_dir <- args

# Validate parameters
if (deconv_method != "xcell") {
  stop("Error: For this analysis, the deconvolution method must be 'xcell'.")
}

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load clinical data
clin_file <- read_tsv(clin_file_path, guess_max = 10000)

# Filter for medulloblastoma subtypes
mb_subtypes <- clin_file %>%
  filter(molecular_subtype %in% c("MB, WNT", "MB, SHH", "MB, Group3", "MB, Group4"))

# Load expression matrix
expr_mat <- readRDS(expr_mat_path)

# Subset expression matrix to include only relevant samples
expr_mat <- expr_mat[, colnames(expr_mat) %in% mb_subtypes]

# Ensure expression matrix has non-zero rows
expr_mat <- expr_mat[rowSums(expr_mat) > 0, ]

# Perform immune cell deconvolution
print("Starting immune cell enrichment analysis using xCell...")
deconv_results <- deconvolute(gene_expression = as.matrix(expr_mat),
                              method = deconv_method,
                              arrays = FALSE)

# Convert results to long format and merge with metadata
deconv_results_long <- deconv_results %>%
  as.data.frame() %>%
  pivot_longer(-cell_type, names_to = "Kids_First_Biospecimen_ID", values_to = "enrichment_score") %>%
  inner_join(mb_subtypes, by = c("Kids_First_Biospecimen_ID"))

# Save results
output_file <- file.path(output_dir, "xcell_immune_enrichment_results.rds")
saveRDS(deconv_results_long, output_file)
print(paste("Enrichment results saved to", output_file))

# Generate Boxplot
immune_enrichment_boxplot <- ggplot(deconv_results_long, aes(x = molecular_subtype, y = enrichment_score, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(), alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "Immune Cell Enrichment by Medulloblastoma Subtypes",
       x = "Molecular Subtype",
       y = "Enrichment Score",
       fill = "Immune Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the boxplot
boxplot_file <- file.path(output_dir, "immune_enrichment_boxplot.png")
ggsave(boxplot_file, immune_enrichment_boxplot, width = 10, height = 6)
print(paste("Boxplot saved to", boxplot_file))

# Generate Sankey Plot
sankey_data <- deconv_results_long %>%
  group_by(molecular_subtype, cell_type) %>%
  summarize(total_enrichment = sum(enrichment_score, na.rm = TRUE)) %>%
  ungroup()

sankey_plot <- ggplot(sankey_data, aes(axis1 = cell_type, axis2 = molecular_subtype, y = total_enrichment)) +
  geom_alluvium(aes(fill = cell_type)) +
  geom_stratum(width = 1/12) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  theme_minimal() +
  labs(title = "Immune Cell Enrichment Sankey Plot",
       x = "Immune Cell Type to Molecular Subtype",
       y = "Enrichment Score",
       fill = "Immune Cell Type")

# Save the Sankey plot
sankey_file <- file.path(output_dir, "immune_enrichment_sankey.png")
ggsave(sankey_file, sankey_plot, width = 12, height = 8)
print(paste("Sankey plot saved to", sankey_file))

# Print completion message
print("Immune cell enrichment analysis and visualization completed!")

