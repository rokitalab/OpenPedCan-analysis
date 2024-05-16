# Gene Set Enrichment Analysis using GSVA

Written by Stephanie J. Spielman to supercede previous analyses in [`ssgsea-hallmark`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/ssgsea-hallmark).
Edited for R 4.4, GSVA 1.52.0, tidyR >1.0 by Jo Lynne Rokita

Primary goals include:

1. Score hallmark pathways based on expression data using GSVA analysis, using a strategy that produces Gaussian-distributed scores.
2. Analyze scores for highly significant differences among tumor classifications 

## Usage:

Note that running this analyis on the full dataset requires > 16GB of memory. 
Run the bash script of this analysis module:

using OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 analyses/gene-set-enrichment-analysis/run-gsea.sh
```

OR by default uses histologies.tsv from data folder
```sh
bash analyses/gene-set-enrichment-analysis/run-gsea.sh
```

*This command above assumes you are in the top directory, OpenPBTA-analysis*

## Folder Content

+ `01-conduct-gsea-analysis.R` performs the GSVA analysis using RSEM TPM expression data for both exome_capture, stranded and polyA data. Results are saved in `results/` TSV files when run via `run-gsea.sh`.

+ `02-model-gsea.Rmd` performs ANOVA and Tukey tests on GSVA scores to evaluate, for each hallmark pathway, differences in GSVA across groups (e.g. short histology or disease type).

+ `results/gsva_scores.tsv` represents GSVA scores calculated from `gene-expression-rsem-tpm-collapsed.rds` with `Rscript --vanilla 01-conduct-gsea-analysis.R`

+ `results/gsva_scores_polya.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` with with `Rscript --vanilla 01-conduct-gsea-analysis.R`

+ Files named as `results/gsva_<tukey/anova>_<all_possible_RNA_library>_<broad_histology/cancer_group)>.tsv` represent results from modeling
	+ Files created with: `Rscript -e "rmarkdown::render('02-model-gsea.Rmd', clean = TRUE, params=list(is_ci = ${IS_CI}))"`
	+ Assumes `results/gsva_scores.tsv` 
 
