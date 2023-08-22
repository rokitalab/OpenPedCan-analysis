
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(analysis_dir, "input")

# Input files
initial_file <- file.path(input_dir, "biospecimen_id_to_bed_map-initial.tsv")
data_file <- file.path(input_dir, "biospecimen_id_to_bed_map.tsv")


initial_df <- readr::read_tsv(initial_file, guess_max = 100000, show_col_types = FALSE)

data_df <- readr::read_tsv(data_file, guess_max = 100000, show_col_types = FALSE)
