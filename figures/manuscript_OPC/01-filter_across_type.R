# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Chante Bethell for CCDL 2019
#

## load libraries
library(tidyverse)
library(ggplot2)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
output_dir <- file.path(root_dir, "figures", "manuscript_OPC")
results_dir <- file.path(output_dir, "results")
plots_dir <- file.path(output_dir, "plots")

# Create directories to hold the output.
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in dataset and remove NAs
histologies_df <-
  readr::read_tsv(file.path(root_dir, "data", "histologies.tsv"), guess_max = 10000) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(harmonized_diagnosis))

# Filter the histologies file to account for multiple samples from the same
# individual and the fact that multiple experimental strategies are in this
# data.frame

# Retain only tumors for this analysis
histologies_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue")

# data.frame with the count of each unique cancer type expression
disease_expression <- histologies_df %>%
  # some recurrences can have different harmonized_diagnosis values
  dplyr::distinct(Kids_First_Participant_ID, harmonized_diagnosis) %>%
  dplyr::group_by(harmonized_diagnosis) %>%
  dplyr::count(name = "count") %>%
  dplyr::arrange(dplyr::desc(count))

# Calculate the total count of the dataset
sum_count <- sum(disease_expression$count)

# Create a percent variable and round to 4 decimal places
# (so values will have 2 decimal places as percentages)
disease_expression <- disease_expression %>%
  dplyr::mutate(percent = paste0((round(count / sum_count, 4) * 100), "%"))

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$harmonized_diagnosis <- with(disease_expression,
                                                reorder(harmonized_diagnosis, -count))

# Write to tsv file
readr::write_tsv(disease_expression,
                 file.path(results_dir,
                           "disease_expression.tsv"))

# Create a bar plot of sample distribution across cancer types
gg_types <- disease_expression %>%
  ggplot2::ggplot(ggplot2::aes(x = harmonized_diagnosis, y = count, fill = count)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Cancer Types", y = "Count",
                title = "Sample Distribution Across Cancer Types") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1200, by = 100)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 75,
    hjust = 1,
    size = 8
  ),
  panel.grid = ggplot2::element_blank()) #+
## percentage is not clear in the figure
  #ggplot2::geom_text(nudge_y = 6.5, size = 2,
                     #ggplot2::aes(label = paste0(disease_expression$percent)))


# Save plot
ggplot2::ggsave(
  gg_types,
  file = file.path(plots_dir, "distribution_across_cancer_types.pdf"),
  width = 22,
  height = 10
)

## get the counts from tumor descriptors
tumor_descriptor <- histologies_df %>%
  dplyr::select(sample_id,
                tumor_descriptor) %>%
  dplyr::distinct() %>%
  dplyr::group_by(tumor_descriptor) %>%
  dplyr::tally() %>%
  dplyr::arrange(dplyr::desc(n))

ggplot(tumor_descriptor, aes(x = factor(tumor_descriptor, levels = tumor_descriptor), y = n, fill = n)) +
  geom_col() +
  theme_bw() +
  labs(x = "Tumor Descriptors", y = "Count", title = "Sample Distribution Across Tumor Descriptors") +
  theme(axis.text.x = element_text(
    angle = 75,
    hjust = 1,
    size = 8
  ),
  panel.grid = element_blank())
ggsave(
  file = file.path(plots_dir, "distribution_across_tumor_descriptors.pdf"),
  width = 22,
  height = 10
)  

