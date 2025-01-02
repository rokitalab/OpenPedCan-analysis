suppressPackageStartupMessages({
  library(optparse)
  library(stats)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
})

# parse parameters
option_list <- list(
  make_option(c("--quantiseq_file"), type = "character",
              help = "immune-deconv results (.rds)"),
  make_option(c("--output_dir"), type = "character", 
              help = "output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))
quan_file <- opt$quantiseq_file
output_dir <- opt$output_dir

# output file
dir.create(output_dir, showWarnings = F, recursive = T)
output_file_df <- file.path(output_dir, paste0( "p_values_output.rds"))
output_file_df1 <- file.path(output_dir, paste0( "p_values_output_accounting_for_cell_type.rds"))
output_plot1 <- file.path(output_dir, "figure1.pdf")

quan <- readRDS(quan_file)

quan_m <- quan %>% 
  filter(cancer_group == "Medulloblastoma")

all <- kruskal.test(fraction ~ molecular_subtype, data = quan_m)
comb <- as.data.frame(combn(unique(quan_m$molecular_subtype), 2))

df <- data.frame(molecular_subtypes=character(),
                 p_value=character(),
                 frd=character(),
                 stringsAsFactors=FALSE)
df[1,1]="all"
df[1,2]=all$p.value
df[1,3]=NA

for (i in c(1:length(comb))){
  df[i+1,1]=paste0(comb[1,i],":",comb[2,i])
  quan_m_sub=quan_m%>%filter(molecular_subtype==comb[1,i] | molecular_subtype==comb[2,i])
  df[i+1,2]=kruskal.test(fraction ~ molecular_subtype, data = quan_m_sub)$p.value
  
}
df$frd <- p.adjust(df$p_value,"BH")

df1 <- data.frame(cell_type=character(),
                  p_value=character(),
                  frd=character(),
                  stringsAsFactors=FALSE)

for(i in seq_along(unique(quan_m_sub$cell_type))){
  df1[i,1]=unique(quan_m_sub$cell_type)[i]
  zz=quan_m_sub %>% filter(cell_type== unique(quan_m_sub$cell_type)[i])
  df1[i,2]=kruskal.test(fraction ~ molecular_subtype, data = zz)$p.value
}


df1$frd <- p.adjust(df1$p_value,"BH")
#select relevant variables
a1 <- quan_m %>% select(Kids_First_Biospecimen_ID,molecular_subtype,fraction, cell_type) %>%
  mutate(name=paste0(Kids_First_Biospecimen_ID,"_",molecular_subtype)) %>% select(-Kids_First_Biospecimen_ID)
  
#convert a1 dataframe for use in PCA function, so that %PC are determined by cell-type counts
a2 <- a1 %>% spread(key="cell_type",
       value= "fraction")

a2 <- as.data.frame(a2)
rownames(a2)<-a2$name
a2 <- a2[,-c(1:2)]
rown <- paste0(sub("^[^ ]* ", "", rownames(a2)),".",c(1:length(rownames(a2))))
rownames(a2)=rown

#PCA analysis
pc <- prcomp(a2)


pca.data <- data.frame(Sample=rownames(pc$x),
                      X=pc$x[,1],
                      Y=pc$x[,2],
                      molecular_subtype=gsub("\\..*","",rownames(pc$x))
                      )
pca.var <- pc$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)

#PCA plot
p1 <- ggplot(data=pca.data, aes(x=X, y=Y, color = molecular_subtype))+
  geom_point()+
  scale_color_manual(values = c("red", "blue", "green", "yellow"))+
  xlab(paste("PC1 _", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 _", pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("Differences in the immune cell distribution between Medulloblastoma subtypes ")+
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(output_plot1, plot = p1, width = 8, height = 9, dpi = 300)
write.csv(df, file = output_file_df)
write.csv(df1, file = output_file_df1)



