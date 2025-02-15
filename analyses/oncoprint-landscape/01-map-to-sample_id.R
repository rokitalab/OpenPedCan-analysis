# J. Taroni for ALSF CCDL 2019
# This script processes MAF, CNV, fusion files and prepares them for oncoprint 
# plotting.
#
# NOTES:
#   * The `Tumor_Sample_Barcode` will now corresponds to the `sample_id` column
#     in the histologies file
#   * We remove ambiguous `sample_id` -- i.e., where there are more than two
#     tumor biospecimens that map to the same sample id
#   * Filtering via an independent specimen file is optional, but highly
#     recommended
#
# EXAMPLE USAGE:
#
# Rscript --vanilla 01-map-to-sample-id.R \
#   --maf_file ../../data/pbta-snv-consensus-mutation.maf.tsv.gz \
#   --cnv_autosomes_file ../../data/consensus_seg_annotated_cn_autosomes.tsv.gz \
#   --cnv_xy_file ../../data/consensus_seg_annotated_cn_x_and_y.tsv.gz \
#   --fusion_file ../../data/pbta-fusion-putative-oncogenic.tsv \
#   --metadata_file ../../data/pbta-histologies.tsv \
#   --palette_file ../../figures/palettes/broad_histology_cancer_group_palette.tsv \
#   --output_directory ../../scratch/oncoprint_files \
#   --filename_lead "primary_only" \
#   --independent_specimens ../../data/independent-specimens.wgswxs.primary.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

#### Command line options ------------------------------------------------------

option_list <- list(
  optparse::make_option(
    c("--maf_file"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains SNV information + hotspots calls",
  ),
  optparse::make_option(
    c("--cnv_autosomes_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains autosome chromosome CNV information"
  ),
  optparse::make_option(
    c("--cnv_xy_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains X and Y chromosome CNV information"
  ),
  optparse::make_option(
    c("--fusion_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information"
  ),
  optparse::make_option(
    c("--metadata_file"),
    type = "character",
    default = NULL,
    help = "file path to clinical file"
  ),
  optparse::make_option(
    c("-p", "--palette_file"),
    type = "character",
    default = NULL,
    help = "file path to the palette file"
  ),  
  optparse::make_option(
    c("--filename_lead"),
    type = "character",
    default = NULL,
    help = "character string used at the beginning of output file names"
  ),
  optparse::make_option(
    c("--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory for files that are prepared in this script"
  ),
  optparse::make_option(
    c("--independent_specimens"),
    type = "character",
    default = NULL,
    help = "optional file path to independent specimens list"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### Directory and output file set up ------------------------------------------

output_dir <- opt$output_directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# All of the output files will use the filename lead supplied as an argument
# followed by what kind of data they contain
maf_output <- file.path(output_dir, paste0(opt$filename_lead, "_maf.tsv"))
fusion_output <- file.path(output_dir, paste0(opt$filename_lead, "_fusions.tsv"))
cnv_output <- file.path(output_dir, paste0(opt$filename_lead, "_cnv.tsv"))

#### Read in data --------------------------------------------------------------
histologies_df <- readr::read_tsv(opt$metadata_file, guess_max = 10000)

palette_df <- readr::read_tsv(opt$palette_file)

# Only keep required maf columns 
keep_cols <-c("Hugo_Symbol", 
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2",
              "Variant_Classification",
              "Variant_Type",
              "Tumor_Sample_Barcode",
              "HGVSp_Short")

maf_df <- data.table::fread(opt$maf_file, skip = 1) %>%
  select(keep_cols)

cnv_autosomes_df <- data.table::fread(opt$cnv_autosomes_file) %>%
  left_join(select(histologies_df, c("Kids_First_Biospecimen_ID","germline_sex_estimate")),
                   by = c("biospecimen_id"="Kids_First_Biospecimen_ID"))
cnv_xy_df <- data.table::fread(opt$cnv_xy_file) 
cnv_df <- rbind(cnv_autosomes_df, cnv_xy_df)
fusion_df <- data.table::fread(opt$fusion_file)

#### Get rid of ambiguous and non-tumor samples --------------------------------

# An ambiguous sample_id will have more than 2 rows associated with it in the
# histologies file when looking at tumor samples -- that means we won't be able
# to determine when an WGS/WXS assay maps to an RNA-seq assay for the purpose of
# the oncoprint plot
ambiguous_sample_ids <- histologies_df %>%
  filter(sample_type == "Tumor",
         composition == "Solid Tissue") %>%
  group_by(sample_id) %>%
  tally() %>%
  filter(n > 2) %>%
  pull(sample_id)

ambiguous_biospecimens <- histologies_df %>%
  filter(sample_id %in% ambiguous_sample_ids) %>%
  pull(Kids_First_Biospecimen_ID)

# we're only going to look at tumor samples in the oncoprint plot
not_tumor_biospecimens <- histologies_df %>%
  filter(sample_type != "Tumor",
         composition != "Solid Tissue") %>%
  pull(Kids_First_Biospecimen_ID)

biospecimens_to_remove <- unique(c(ambiguous_biospecimens,
                                   not_tumor_biospecimens))

# Filter the files!
maf_df <- maf_df %>%
  dplyr::filter(!(Tumor_Sample_Barcode %in% biospecimens_to_remove))
cnv_df <- cnv_df %>%
  dplyr::filter(!(biospecimen_id %in% biospecimens_to_remove))
fusion_df <- fusion_df %>%
  dplyr::filter(!(Sample %in% biospecimens_to_remove))
histologies_df <- histologies_df %>%
  dplyr::filter(!(Kids_First_Biospecimen_ID %in% biospecimens_to_remove))

#### Filter to independent specimens (optional) --------------------------------

if (!is.null(opt$independent_specimens)) {

  ind_biospecimen <- readr::read_tsv(opt$independent_specimens) %>%
    pull(Kids_First_Biospecimen_ID)

  # filter the genome data, e.g., SNV and CNV data, to only include biospecimen
  # identifiers in the independent file
  maf_df <- maf_df %>%
    filter(Tumor_Sample_Barcode %in% ind_biospecimen)
  cnv_df <- cnv_df %>%
    filter(biospecimen_id %in% ind_biospecimen)

  # for the RNA-seq samples, we need to map from the sample identifier
  # associated with the independent specimen and back to a biospecimen ID
  ind_sample_id <- histologies_df %>%
    filter(Kids_First_Biospecimen_ID %in% ind_biospecimen) %>%
    pull(sample_id)

  # get the corresponding biospecimen ID to be used
  rnaseq_ind <- histologies_df %>%
    filter(sample_id %in% ind_sample_id,
           experimental_strategy == "RNA-Seq") %>%
    pull(Kids_First_Biospecimen_ID)

  # finally filter the fusions
  fusion_df <- fusion_df %>%
    filter(Sample %in% rnaseq_ind)

  # filter histologies file
  histologies_df <- histologies_df %>%
    filter(Kids_First_Biospecimen_ID %in% c(ind_biospecimen,rnaseq_ind))

}

#### MAF file preparation ------------------------------------------------------

message("Preparing MAF file...")

# join the sample_id information to the MAF file and then set as the tumor
# sample barcode
maf_df <- maf_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  # now let's remove this `Tumor_Sample_Barcode` column with biospecimen IDs in
  # preparation for our next step -- renaming `sample_id`
  select(-Tumor_Sample_Barcode) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Write MAF to file
readr::write_tsv(maf_df, maf_output)

#### Fusion file preparation ---------------------------------------------------

message("Preparing fusion file...")

# We'll handle fusions where reciprocal fusions exist (e.g., 
# reciprocal_exists == TRUE) separately from other fusions
fusion_reciprocal_df <- fusion_df %>%
  # Reciprocal fusions only
  filter(reciprocal_exists) %>%
  # BSID + Gene1--Gene2
  select(Sample, FusionName) %>%
  # Because we're only looking at presence or absence here, we can filter to 
  # distinct identifier-fusion pairs
  distinct() %>%
  group_by(Sample, FusionName) %>%
  # Put genes in the fusions in alphabetical order to collapse
  mutate(SortedFusionName = str_c(
    sort(str_split(FusionName,
                   "--",
                   simplify = TRUE),
    ), 
    collapse = "--")
  ) %>%
  ungroup() %>%
  select(Sample,
         FusionName = SortedFusionName) %>%
  # When fusions are in alphabetical order, this ensures each reciprocal
  # fusion is only counted once
  distinct()

# No need to put fusions where no reciprocal exists in alphabetical order,
# so we handle these separately
fusion_no_reciprocal_df <- fusion_df %>%
  filter(!reciprocal_exists) %>%
  select(Sample, FusionName) %>%
  # But we can remove duplicates to avoid them being counted as multihit
  distinct()

# Bind reciprocal and no reciprocal together
fusion_filtered_df <- bind_rows(fusion_reciprocal_df, fusion_no_reciprocal_df)
rm(fusion_reciprocal_df, fusion_no_reciprocal_df)

# Separate fusion gene partners
fus_sep <- fusion_filtered_df %>%
  # Separate the 5' and 3' genes
  tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--") %>%
  # Use row numbers to mark unique fusions - this will help us when
  # we melt and remove selfie fusions below
  tibble::rowid_to_column("Fusion_ID") %>%
  select(Fusion_ID, Sample, Gene1, Gene2) %>%
  reshape2::melt(id.vars = c("Fusion_ID", "Sample"),
                 variable.name = "Partner",
                 value.name = "Hugo_Symbol") %>%
  arrange(Fusion_ID)

multihit_fusions <- fus_sep %>%
  # For looking at multi-hit fusions, we want to remove selfie fusions
  # If we drop the 5', 3' information in Partner, we can use distinct
  # to remove the selfie fusions because Fusion_ID marks what fusion
  # it came from
  select(-Partner) %>%
  distinct() %>%
  # Now we want to count how many times a gene is fused within a sample
  # Anything with more than one count is considered multi-hit
  group_by(Sample, Hugo_Symbol) %>%
  dplyr::count() %>%
  filter(n > 1) %>%
  select(-n) %>%
  mutate(Variant_Classification = "Multi_Hit_Fusion")

# Filter out multi-hit fusions from the other fusions that we will
# label based on whether they are the 5' or 3' gene
single_fusion <- fus_sep %>%
  anti_join(multihit_fusions,
            by = c("Sample", "Hugo_Symbol")) %>%
  select(Sample, Hugo_Symbol) %>%
  distinct() %>%
  mutate(Variant_Classification = "Fusion")

# Combine multi-hit and single fusions into one data frame
# And add the other identifiers!
reformat_fusion <- multihit_fusions %>%
  bind_rows(single_fusion) %>%
  mutate(Variant_Type = "OTHER") %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id,
         Kids_First_Biospecimen_ID = Sample)

# Write to file
readr::write_tsv(reformat_fusion, fusion_output)

#### CNV file preparation ------------------------------------------------------

message("Preparing CNV file...")
cnv_df <- cnv_df %>%
  inner_join(select(histologies_df,
                    Kids_First_Biospecimen_ID,
                    sample_id),
             by = c("biospecimen_id"="Kids_First_Biospecimen_ID")) %>%
  mutate(Tumor_Sample_Barcode =  sample_id) %>%
  dplyr::rename(Variant_Classification = status,
         Hugo_Symbol = gene_symbol) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
  # mutate loss and amplification to Del and Amp to fit Maftools format
  dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "deep deletion" ~ "Del",
                                                          Variant_Classification == "amplification" ~ "Amp",
                                                          TRUE ~ as.character(Variant_Classification))) %>%
  # only keep Del and Amp calls
  filter(Variant_Classification %in% c("Del", "Amp"))

# Write to file
readr::write_tsv(cnv_df, cnv_output)

message("Gathering N per cancer_group...")
# Gather Ns
gather_n_per_cancer_group <- function(metadata, bh){
  if (bh != "Other CNS") {
    metadata <- metadata %>%
      dplyr::filter(broad_histology == bh) %>%
      select(cancer_group_display, sample_id) %>%
      unique() %>%
      group_by(cancer_group_display) %>%
      tally() %>%
      as.data.frame()

  } else {
    metadata <- metadata %>%
      dplyr::filter(
        broad_histology %in% c(
          "Ependymal tumor",
          "Tumors of sellar region",
          "Neuronal and mixed neuronal-glial tumor",
          "Tumor of cranial and paraspinal nerves",
          "Meningioma",
          "Mesenchymal non-meningothelial tumor",
          "Germ cell tumor",
          "Choroid plexus tumor",
          "Histiocytic tumor",
          "Tumor of pineal region",
          "Metastatic tumors",
          "Other astrocytic tumor",
          "Lymphoma",
          "Melanocytic tumor",
          "Other tumor"
        )
      )  %>%
      select(cancer_group_display, sample_id) %>%
      unique() %>%
      group_by(cancer_group_display) %>%
      tally() %>%
      as.data.frame()

  }
}

# Add cancer_group_display into histologies_df
histologies_df <- histologies_df %>%
  inner_join(
    select(palette_df, broad_histology, cancer_group, cancer_group_display)
  ) %>%
  # Remove rows with NA in cancer_group_display
  tidyr::drop_na(cancer_group_display)

output_file <- file.path('tables', paste0(opt$filename_lead, "_n_per_cancer_group.tsv"))
lapply(as.list(c("Low-grade astrocytic tumor","Embryonal tumor",
         "Diffuse astrocytic and oligodendroglial tumor",
         "Other CNS")), function(x) gather_n_per_cancer_group(histologies_df,x)) %>%
	bind_rows() %>%
        unique() %>%
        readr::write_tsv(output_file)

