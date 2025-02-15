# Author: Sangeeta Shukla
# This script servers as a precursor to the DESeq analysis step, as it calculates the GTEx_index, Hist_Index values
# This script also creates Histology and Counts data subsets which satisfy given clinical criteria 


# Load required libraries

suppressPackageStartupMessages({
  library(optparse)
  library(sys)
})


option_list <- list(
    make_option(c("-c", "--hist_file"), type = "character",
              help = "Histology data file (.TSV)"),
    make_option(c("-n", "--counts_file"), type = "character",
              help = "Gene Counts file (.rds)"),
    make_option(c("-o", "--outdir"), type = "character",
              help = "Path to Output Directory", default = "."),
    make_option(c("-y", "--ind_allcohorts"), type = "character",
              help = "Independent specimens list for all cohorts(.tsv)"),
    make_option(c("-z", "--ind_eachcohort"), type = "character",
              help = "Independent specimens list for each cohort (.tsv)")
)


opt <- parse_args(OptionParser(option_list = option_list))


#Load histology file
hist <- read.delim(opt$hist_file, header=TRUE, sep = '\t')

#Load expression counts data
countData <- readRDS(opt$counts_file)

#Load list of independent specimens for across all cohorts
ind_spec_all_cohorts <- read.delim(opt$ind_allcohorts, header=TRUE, sep='\t')

#Load list of independent specimens for each cohort
ind_spec_each_cohort <- read.delim(opt$ind_eachcohort, header=TRUE, sep='\t')

outdir <- opt$outdir
cmd_mkdir <- paste("mkdir",outdir,sep=" ")
system(cmd_mkdir)


# Subset Histology file for samples only found in the current the countData file (To ensure no discepancies cause errors later in the code)
hist.filtered <- unique(hist[which(hist$Kids_First_Biospecimen_ID %in%  colnames(countData)),])


# Create an array of unique cancer_group found in histologies.tsv
cancerGroup <- unique(hist.filtered$cancer_group)
cancerGroup <- cancerGroup[which(!is.na(cancerGroup))]

# Create an array of unique research cohorts found in histologies.tsv
resCohort <- unique(hist.filtered$cohort)
resCohort <- resCohort[which(!is.na(resCohort))]

# Combine the cancer_group and cohort as columns in a new array
cancerGroup_cohort_set <- expand.grid(cancerGroup=cancerGroup,cohort=resCohort)

# Create a new array which can take each combination of cancer_group+cohort 
# Add another column with counts of patients whose data is available for that combination
patientCount_set <- data.frame()

for (I in 1:length(cancerGroup_cohort_set$cancerGroup))
{
  patientCount_set <- rbind(patientCount_set, 
                            data.frame(cancerGroup=cancerGroup_cohort_set$cancerGroup[I], 
                                       cohort=cancerGroup_cohort_set$cohort[I], 
                                       counts=length(unique(hist.filtered$Kids_First_Biospecimen_ID[
                                         which(hist.filtered$cancer_group == cancerGroup_cohort_set$cancerGroup[I] 
                                               & hist.filtered$cohort == cancerGroup_cohort_set$cohort[I])
                                       ]
                                       )
                                       )
                            )
  )
}


patientCount_set <- subset(patientCount_set,patientCount_set$counts>=3)

hist.filtered_final = data.frame()
for(K in 1:nrow(patientCount_set))
{
  hist.filtered_final <- rbind(hist.filtered_final,hist.filtered[which(hist.filtered$cancer_group == patientCount_set$cancerGroup[K] &
                                                                          hist.filtered$cohort == patientCount_set$cohort[K] ),])
}
hist.filtered_final <- rbind(hist.filtered_final,hist.filtered[which(!is.na(hist.filtered$gtex_group)),])





#Subset countdata for data that are present in the hitstoly files (To ensure no discepancies cause errors later in the code)
countData_filtered <- countData[,which(colnames(countData) %in% hist.filtered_final$Kids_First_Biospecimen_ID)]


#Save all the unique cancer histologies in a variable. These cancer histologies represent the patient data in the countsdata
Cancer_Histology <- unique(hist.filtered_final$cancer_group)

#Save all the GTEx tissue subgroups in a variable. These cancer histologies represent the GTEx RNDA data available in the countsdata
Gtex_Tissue_subgroup <- sort(unique(hist.filtered_final$gtex_subgroup))


#Save all the cohorts represented in the countsdata into a variable. Remove all 'NA's from the list. 
#And paste cohort to cancer group (eg GMKF_Neuroblastoma)
Cancer_Histology_COHORT <- unique(
  paste(hist.filtered_final$cohort[which(!is.na(hist.filtered_final$cancer_group))],
        hist.filtered_final$cancer_group[which(!is.na(hist.filtered_final$cancer_group))],
        sep="_")
)

#Save all the histologies represented in the countsdata into a variable. 
#Remove all 'NA's from the list. 
#This will be the basis of all the data from each histology combined regardless of cohort (eg all-cohorts_Neuroblastoma)
Cancer_Histology <- paste("all-cohorts",Cancer_Histology[which(!is.na(Cancer_Histology))],sep="_")


#Save all the GTEx subgroups represented in the countsdata into a variable. Remove all 'NA's 
Gtex_Tissue_subgroup <- Gtex_Tissue_subgroup[!is.na(Gtex_Tissue_subgroup)]

#Create an empty df to populate with rbind of all normal Kids_First_Biospecimen_ID and gtex_subgroup
#Create DF that list all Kids_First_Biospecimen_IDs by GTEX subgroup
sample_type_df_normal <- data.frame()
for(I in 1:length(Gtex_Tissue_subgroup))
{
  sample_type_df_normal <- rbind(sample_type_df_normal,
                                 data.frame(Case_ID = hist.filtered_final$Kids_First_Biospecimen_ID[which(hist.filtered_final$gtex_subgroup == Gtex_Tissue_subgroup[I])]
                                            ,Type = Gtex_Tissue_subgroup[I]), stringsAsFactors = FALSE)
}

#Create an empty df to populate with rbind of all tumor Kids_First_Biospecimen_ID and cancer_group
#Create DF that list all Kids_First_Biospecimen_IDs by cancer group subgroup
sample_type_df_tumor <- data.frame()
for(I in 1:length(Cancer_Histology))
{
  sample_type_df_tumor <- rbind(sample_type_df_tumor,data.frame(Case_ID = hist.filtered_final$Kids_First_Biospecimen_ID[which(hist.filtered_final$cancer_group == gsub("all-cohorts_","",Cancer_Histology[I]))]
                                                                ,Type=Cancer_Histology[I], stringsAsFactors = FALSE))
}

#Only use samples that are independent/unique across all cohorts
sample_type_df_tumor <- subset(sample_type_df_tumor,sample_type_df_tumor$Case_ID %in% ind_spec_all_cohorts$Kids_First_Biospecimen_ID)



#Create an empty df to populate with rbind of all tumor Kids_First_Biospecimen_ID and cancer_group by cohort
#Create DF that list all Kids_First_Biospecimen_IDs by Cohort - Cancer groups
sample_type_df_tumor_cohort <- data.frame()
for(I in 1:length(Cancer_Histology_COHORT))
{
  Cancer_Histology_COHORT_cohort <- strsplit(Cancer_Histology_COHORT[I],split="_")[[1]][1]
  Cancer_Histology_COHORT_cancer_group <- strsplit(Cancer_Histology_COHORT[I],split="_")[[1]][2]
  sample_type_df_tumor_cohort <- rbind(sample_type_df_tumor_cohort,
                                       data.frame(Case_ID = hist.filtered$Kids_First_Biospecimen_ID[which(hist.filtered$cancer_group == Cancer_Histology_COHORT_cancer_group 
                                                                                                          & hist.filtered$cohort == Cancer_Histology_COHORT_cohort)]
                                                  ,Type=Cancer_Histology_COHORT[I], stringsAsFactors = FALSE))
}


#Only use samples that are independent/unique for each cohort
sample_type_df_tumor_cohort <- subset(sample_type_df_tumor_cohort,sample_type_df_tumor_cohort$Case_ID %in% ind_spec_each_cohort$Kids_First_Biospecimen_ID)


#Combine the rows from the normal and tumor sample df
sample_type_df <- rbind(sample_type_df_tumor,sample_type_df_tumor_cohort,sample_type_df_normal)


#Filter one more to ensure the rownames in the countsdata file match the sample dataframe for DEG just created
countData_filtered_DEG <- countData_filtered[,which(colnames(countData_filtered) %in% sample_type_df$Case_ID)]

#Filter one more to ensure the rownames in the countsdata file match the sample dataframe for DEG just created
sample_type_df_filtered <- unique(sample_type_df[which(sample_type_df$Case_ID %in% colnames(countData_filtered_DEG)),])

#Define All cancer groups (Combined and cohort-specific) in the histology list
histology_filtered <- unique(sample_type_df_filtered$Type[-grep("^GTEX",sample_type_df_filtered$Case_ID)])

#Define All GTEx groups as normal in the GTEX_filtered list
GTEX_filtered <- unique(sample_type_df_filtered$Type[grep("^GTEX",sample_type_df_filtered$Case_ID)])



#Define all expected filenames as pairs of histology tissue and cancer_group
filename_list <- data.frame()
for(I in 1:length(histology_filtered))
{
  for (J in 1:length(GTEX_filtered))
  {
    FILENAME <- gsub("all-cohorts","all_cohorts",gsub(" |/|;|:|\\(|)","_",paste(histology_filtered[I],GTEX_filtered[J],sep="_v_")))
    filename_list <- rbind(filename_list, FILENAME)
  }
}



# Print output to files
fileConn_GTEx<-file(paste(outdir,"/GTEx_Index_limit.txt",sep=""),open = "w")
write.table(length(GTEX_filtered), file = fileConn_GTEx, append = FALSE, row.names = FALSE, col.names = FALSE)
close(fileConn_GTEx) 

fileConn_Hist<-file(paste(outdir,"/Hist_Index_limit.txt",sep=""),open = "w")
write.table(length(histology_filtered), file = fileConn_Hist, append = FALSE, row.names = FALSE, col.names = FALSE)
close(fileConn_Hist) 


fileConn_FilenameList<-file(paste(outdir,"/filelist_expected.txt",sep=""),open = "w")
write.table(filename_list, file = fileConn_FilenameList, append = FALSE, row.names = FALSE, col.names = FALSE)
close(fileConn_FilenameList)


write.table(hist.filtered_final, file=paste(outdir,"/histologies_subset.tsv",sep=""), sep="\t", col.names = T, row.names = F,quote = F)
saveRDS(countData_filtered,file=paste(outdir,"/countData_subset.rds",sep=""))


Hist_GTEX_Indices = expand.grid(1:length(GTEX_filtered),1:length(histology_filtered))[,c(2,1)]
write.table(Hist_GTEX_Indices, file = "indices.txt", sep = "\t",quote = F, row.names = F, col.names = F)





