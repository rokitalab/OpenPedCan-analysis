# Script will take the outputs from 01-immune-deconv.R and investigate immune differences between different types of medullobstastoma subtypes.
# This script will use the xcell outputs since they are the preferred deconvolution method.
# Author: Audrey Crowther

# Load Packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the xcell input file
xcell_data <- readRDS("~/OpenPedCan-Project/analyses/immune-deconv/results/xcell_output.rds")

# Subset to get only the Medullablastoma tumors
medullablastoma <- xcell_data %>%
  filter(cancer_group == "Medulloblastoma")

# Make subtypes names cleaner. Shortened names will be used in data visualization plots
medullablastoma$molecular_subtype <- gsub("MB, ", "", medullablastoma$molecular_subtype)

# Ensure all values are read as numeric
medullablastoma$fraction <- as.numeric(medullablastoma$fraction)

# A MANOVA (Multivariate ANOVA) will be used to assess immune composition differences - https://www.appsilon.com/post/manova-in-r
# While an anova will investigate if there are immune cell differences between subtypes on a per immune cell basis,
# there are many immune cell types in the data file, and the fraction of each immune cell is presumably dependent on one another 
# (e.g. fractions are presumably out of 100 for all tissue types in the tumor).
# There is an assumption of normality with the MANOVA, similar to an ANOVA, which may not be appropriate for a given dataset.
# Alternatively, FDR-correcting for repeated ANOVAs may also suffice. No FDR correction may be appropriate for most exploratory option for immune differences in medullablastoma tumors.

# Reshape so that the immune cells are columns. 
# Also only keeps the patient IDs and molecular subtype column for simplicity.
medullablastoma_reshape<-pivot_wider(medullablastoma,
                                  id_cols = c("Kids_First_Biospecimen_ID", "molecular_subtype"),
                                  names_from = "cell_type",
                                  values_from = "fraction",
                                  values_fill = FALSE)

# Isolate the immune cell fractions 
immune_cell_types <- medullablastoma_reshape[, 3:ncol(medullablastoma_reshape)]


# Fit MANOVA model
manova_result <- manova(as.matrix(immune_cell_types) ~ molecular_subtype, data = medullablastoma_reshape)

# summary(manova_result) is throwing an error. There are no detected NA's in the input matrix.
# Potentially caused by correlated variables - which would make sense for the data structure 
# e.g. T cell fraction will likely correlate with CD4+ fraction and further subsets. 
# Adding further requirements for filtering/data processing can be added later.  


# Function to perform post-hoc analysis for manova results. 
# This will specify immune differences per medullablastoma subtype in a pair-wise comparison. 
extract_posthoc <- function(var) {
  # Enclose column name in backticks to handle spaces
  result <- TukeyHSD(aov(as.formula(paste0("`", var, "` ~ molecular_subtype")), data = medullablastoma_reshape))
  data.frame(result$molecular_subtype) %>%
    mutate(immune_cell_type = var, comparison = rownames(result$molecular_subtype)) %>%
    select(immune_cell_type, comparison, everything())
}

# Apply post-hoc analysis to each immune cell type
immune_cell_names <- colnames(immune_cell_types)

# Combine results into one matrix
posthoc_results_combined <- lapply(immune_cell_names, extract_posthoc)

# Combine all post-hoc results into a single data frame
posthoc_results_df <- bind_rows(posthoc_results_combined)

# Sort results by p-value (most significant first)
sorted_results <- posthoc_results_df %>% arrange(p.adj)

output_path <- "~/OpenPedCan-Project/analyses/immune-deconv/results/Medullablastoma_ImmuneEnrichments.rds"

# Save results file
saveRDS(sorted_results, output_path)


# Function to create boxplots for each immune cell type
plot_immune_cell <- function(cell_type, data) {
  plot <- ggplot(data, aes(x = molecular_subtype, y = .data[[cell_type]], fill = molecular_subtype)) +
    labs(y = "Immune Cell Fraction", x = NULL) +
    ggtitle(cell_type) +
    geom_boxplot(outlier.shape = NA) +  # Remove outliers
    geom_jitter(shape = 1, position = position_jitter(0.2), size = 3) +  # Add jitter for individual points
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 13),  # Increase x-axis tick mark labels
      axis.text.y = element_text(size = 13),  # Increase y-axis tick mark labels
      axis.title.x = element_text(size = 13, margin = margin(t = 10, unit = "pt")),  # Increase x-axis title size
      axis.title.y = element_text(size = 14, margin = margin(r = 14, unit = "pt")),  # Increase y-axis title size
      plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 20, unit = "pt"))  # Center and style title
    )
  
  # filename <- gsub( " ", "_", cell_type) # Optional to implement. Change "cell_type" below to "filename"
  
  # Save the plot as a PDF
  ggsave(filename = paste0(cell_type, ".pdf"), plot = plot, 
         path = "~/OpenPedCan-Project/analyses/immune-deconv/results", 
         width = 8, height = 6, units = "in")
}

# Loop through immune cell types and generate plots
immune_cell_names <- colnames(immune_cell_types)
plots <- list()

for (cell_type in immune_cell_names) {
  plot_immune_cell(cell_type, medullablastoma_reshape)
}