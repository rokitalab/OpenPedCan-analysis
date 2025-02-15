#This script assigns ATRT into three known subtypes using methylation result.
#Subtypiong esults is saved as ATRT-molecular-subtypes.tsv

# Set up library
library(tidyverse)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to result
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "results")

# Directory of json file for subseting
atrt_subset_file <- file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "atrt-subset", "ATRT_subtyping_path_dx_strings.json")

# Read in histologies
histo <- 
  readr::read_tsv(file.path(root_dir, "data", "histologies-base.tsv"), guess_max = 100000)

# Read json file
path_dx_list_atrt <- jsonlite::fromJSON(atrt_subset_file)
  
# Filter histo, 
# select all ATRT biospecimens from PBTA and/DGD
atrt_df <- histo %>%
  dplyr::filter(cohort %in% c("PBTA", "Kentucky", "DGD", "PPTC")) %>%
  dplyr::filter(pathology_diagnosis %in% path_dx_list_atrt$exact_path_dx)

# Create a dataframe with sample_id and matched biospecimens id for all expt strategies
atrt_df_methyl <- atrt_df %>% 
  # pull either high-conf DKFZ subtypes OR high-conf NIH subtypes
  filter(experimental_strategy == "Methylation" & 
           (dkfz_v12_methylation_subclass_score >= 0.8 & 
              grepl("ATRT_", dkfz_v12_methylation_subclass)) |
          (NIH_v2_methylation_Superfamily_mean_score >= 0.9 & 
             NIH_v2_methylation_Class_mean_score >= 0.9 & 
             grepl("ATRT", NIH_v2_methylation_Class))
          ) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, match_id, sample_id, composition, 
         dkfz_v12_methylation_subclass, dkfz_v12_methylation_subclass_score,
         NIH_v2_methylation_Class, NIH_v2_methylation_Superfamily_mean_score, NIH_v2_methylation_Class_mean_score) %>%
  # mutate ATRT subtypes from methyl so they are in the format needed later
  mutate(molecular_subtype = case_when(dkfz_v12_methylation_subclass_score >= 0.8 ~ dkfz_v12_methylation_subclass,
                                                   TRUE ~ NA_character_)) %>%
  # now if any are NA, then take the NIH classification
  mutate(molecular_subtype = case_when(is.na(molecular_subtype) & 
                                         NIH_v2_methylation_Superfamily_mean_score >= 0.9 & 
                                         NIH_v2_methylation_Class_mean_score >= 0.9 ~ NIH_v2_methylation_Class,
                                       TRUE ~ molecular_subtype)) %>%
  # rename
  mutate(molecular_subtype = gsub("^ATRT_", "ATRT, ", molecular_subtype))

methyl_no_subtype <- atrt_df %>% 
  filter(experimental_strategy == "Methylation",
         !match_id %in% atrt_df_methyl$match_id) %>% 
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, sample_id, composition, 
         dkfz_v12_methylation_subclass, dkfz_v12_methylation_subclass_score) %>%
  mutate(molecular_subtype = "ATRT, To be classified") 


# make a methyl map
methyl_map <- atrt_df_methyl %>%
  bind_rows(methyl_no_subtype) %>%
  select(match_id, molecular_subtype) %>%
  unique() %>%
  mutate(molecular_subtype_methyl = molecular_subtype)

# any dups? no
length(unique(methyl_map$match_id)) == length(methyl_map$match_id)


# for the samples, whose dkfz_v12_methylation_subclass_score >= 0.8 and dkfz_v12_methylation_subclass is one of the three types in ATRT_subtype_list, 
# their molecular subtypes are same as dkfz_v12_methylation_subclass
# for the samples without methylation sequencing, their molecular subtype are "ATRT, To be classified."
# For the samples fit all the other situations, their molecular subtype are "ATRT, To be classified."
atrt_subtype_final <- atrt_df %>% 
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, match_id, sample_id, composition) %>%
  left_join(methyl_map) %>%
  mutate(molecular_subtype = case_when(is.na(molecular_subtype) ~ "ATRT, To be classified",
                                       TRUE ~ molecular_subtype)) %>%
  # write result
  readr::write_tsv(file.path(results_dir, "ATRT-molecular-subtypes.tsv"))


# how many tumors subtyped?
tumors_subtyped <- atrt_subtype_final %>%
  select(match_id, molecular_subtype) %>%
  unique() %>%
  count(molecular_subtype)

print(tumors_subtyped)
