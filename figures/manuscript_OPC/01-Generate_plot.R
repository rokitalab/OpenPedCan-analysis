# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Chante Bethell for CCDL 2019
# Zhuangzhuang Geng, D3B 2024

## load libraries
library(tidyverse)
library(ggplot2)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
analysis_dir <- file.path(root_dir, "figures", "manuscript_OPC")
plots_dir <- file.path(analysis_dir, "plots")

## source theme_Publication
source(file.path(analysis_dir, "utils", "theme_for_plot.R"))

# Create directories to hold the output.
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in dataset and remove NAs
histologies_df <-
  readr::read_tsv(file.path(root_dir, "data", "histologies.tsv"), guess_max = 10000) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(harmonized_diagnosis))

## Generate the pie plot for tumor description
tumor_descriptor <- histologies_df %>% 
  select(match_id, tumor_descriptor) %>% 
  distinct() %>% 
  group_by(tumor_descriptor) %>% 
  tally() %>% 
  filter(!(tumor_descriptor %in% c("Unknown", NA))) %>% 
  arrange(desc(n))

ggplot(tumor_descriptor, aes(x = "", y = n, fill = factor(tumor_descriptor, levels = tumor_descriptor))) + 
  geom_bar(stat = "identity", width = 1, color = "white") + 
  coord_polar("y", start = 0) + 
  theme_void() + 
  theme(legend.title = element_blank())

ggsave(file.path(plots_dir, "tumor_descriptor.png"), height = 8, width = 10)

## generate stacked barplot for experimental strategy
exp_strategy_df <- histologies_df %>% 
  mutate(group = case_when(cohort %in% c("TCGA", "GTEx") ~ "Adult", TRUE ~ "Pediatric")) %>%
  mutate(group = factor(group, levels = c("Pediatric", "Adult"))) %>%
  select(group, cohort, experimental_strategy)

exp_level <- unique(exp_strategy_df$experimental_strategy)[order(unique(exp_strategy_df$experimental_strategy))]
    
ggplot(exp_strategy_df, aes(x = cohort, fill = factor(experimental_strategy, levels = exp_level))) + 
  geom_bar() +
  scale_fill_manual("Experimental Strategy", values = c("darkgoldenrod2", "yellow1", "magenta2", "salmon", "darkslategray2", "blueviolet", "deeppink", "darkseagreen")) +
  facet_wrap(vars(group)) + 
  theme_Publication()

ggsave(file.path(plots_dir, "experimental_strategy.png"), height = 10, width = 16)

