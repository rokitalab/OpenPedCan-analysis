# In this script we will be gathering pathology diagnosis
# and pathology free text diagnosis terms to select HGG
# samples for downstream HGG subtyping analysis and save 
# the json file in hgg-subset folder

library(tidyverse)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_file <- file.path(root_dir,
                         "analyses",
                         "molecular-subtyping-HGG",
                         "hgg-subset",
                         "hgg_subtyping_path_dx_strings.json")

# Read histologies_base.tsv file
histo <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv")) %>% 
  dplyr::filter(cohort %in% c("PBTA", "Kentucky", "DGD", "PPTC"))

# The `pathology_diagnosis` fields for HGG
# as we identified in 00-v9-HGG-select-pathology-dx.Rmd are:
exact_path_dx<- c(
  "High-grade glioma/astrocytoma (WHO grade III/IV)",
  "Brainstem glioma- Diffuse intrinsic pontine glioma",
  "Glioblastoma",
  "Astrocytoma;Oligoastrocytoma",
  "Astrocytoma",
  "High-grade glioma/astrocytoma (WHO grade III/IV);Oligodendroglioma",
  "High-grade glioma/astrocytoma (WHO grade III/IV);Neurofibroma/Plexiform",
  "Gliomatosis Cerebri;High-grade glioma/astrocytoma (WHO grade III/IV)"
)

# Gliomatosis Cerebri can be high grade glioma or low grade 
# glioma so we will add an inclusion criteria for v9 release 
# to only keep `Gliomatosis Cerebri` samples if pathology_free_text_diagnosis
# as `anaplastic gliomatosis cerebri (who grade 4)`
path_free_text_exact <- c("anaplastic gliomatosis cerebri (who grade 4)",
                                      "astroblastoma")

#Identify which samples are IHGs and bring them into the HGG module for subtyping. 
#This can be done by either searching 
#pathology_free_text_diagnosis for the term/terms: "infant type hemispheric glioma" or 

IHG_path_free_path_dx <- histo %>%
  filter(grepl("infant type hemispheric glioma", pathology_free_text_diagnosis)) %>%
  pull(pathology_free_text_diagnosis) %>%
  unique() 

# create an exclude list - LCH has BRAF V600E, so may come in by mistake
exclude_dx <- c("Langerhans Cell histiocytosis")

# Create a list with the strings we'll use for inclusion.
terms_list <- list(exact_path_dx = exact_path_dx,
                   path_free_text_exact = path_free_text_exact, 
                   IHG_path_free_path_dx = IHG_path_free_path_dx,
                   exclude_path_dx = exclude_dx)


#Save this list as JSON.
writeLines(jsonlite::prettify(jsonlite::toJSON(terms_list)), output_file)
