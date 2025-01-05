## Immune Deconvolution

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi)), updated Kelsey Keith ([@kelseykeith](https://github.com/kelseykeith))

### Description

The goal of this analysis is to use the R package `immunedeconv` to quantify and compare various immune cell types in the tumor microenvironment (TME) across various cancer and GTEx groups. 
The package `immunedeconv`, provides the following deconvolution methods: `xCell` (n = 64; immune and non-immune cell types), `CIBERSORT` (relative mode; n = 22 immune cell types); `CIBERSORT (abs.)` (absolute mode; n = 22 immune cell types), `TIMER` (n = 6), `EPIC` (n = 6), `quanTIseq` (n = 10) and `MCP-Counter` (n = 8). 

Both `CIBERSORT` and `CIBERSORT (abs.)` require two files i.e. `LM22.txt` and `CIBERSORT.R`, that are available upon request from https://cibersort.stanford.edu/. Please refer to https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html#special-case-cibersort for more details. We recommend using `xCell` instead when these files are not available to the user. 

### Method selection


We use two methods: xCell and quanTIseq. 


We chose xCell because it: 
1) is the most comprehensive deconvolution method and is able to deconvolute the maximum number of immune and non-immune cell types 
2) is highly robust against background predictions and 
3) can reliably identify the presence of immune cells at low abundances (0-1% infiltration depending on the immune cell type).

xCell outputs immune scores as arbitrary scores that represent cell type abundance. 
Importantly, these scores may be compared between samples (inter-sample comparisons), but _may not_ be compared across cell types or cancer types, as described in the [`immunedeconv` documentation](https://omnideconv.org/immunedeconv/articles/immunedeconv.html#interpretation-of-scores). This is in part because xCell is actually a signature-based method and not a deconvolution method, as is described in the [xCell Publication](https://doi.org/10.1186/s13059-017-1349-1):
> Unlike signature-based methods, which output independent enrichment scores per cell type, the output from deconvolution-based methods is the inferred proportions of the cell types in the mixture.

Therefore, we also use `quanTIseq` as a complementary method. Although `quanTIseq` looks at fewer cell types, the scores can be interpreted as absolute fractions, thereby allowing comparison _both_ across samples and cell types, [as described](https://omnideconv.org/immunedeconv/articles/immunedeconv.html#interpretation-of-scores).


### Custom Analysis Scripts - 2 Seperate analysis were performed on each deconv_method!

#### AKB - 2025-01 ####

#### 02-immune-deconv.R

1. Inputs from data download by running the following script

```
bash scripts/download-testing-files.sh

data/testing/gene-expression-rsem-tpm-collapsed.rds
data/testing/histologies.tsv
```

```
### Script Changes:

#### Script : 02-Immune-deconv.R

1. An opt parser (cancer_type_filter) has been added to work on any cancer group. In this case, we are filtering for "Medulloblastoma"
2. Plotting for Median Cell fractions and Total Cell fractions for the Cancer molecular subtypes (After filtering, "MB, Group4" and "MB, SHH" were only molecular subtypes found!)

```

2. Running the analysis:

```
bash run-02-immune-deconv.sh
```


3. Output: 

```
For (deconv_methods: xcell and quantiseq)

results_custom_deconv/{deconv_method}_output.rds
results_custom_deconv/{deconv_method}_output_cellfractions.pdf
```


#### 03-Immune-deconv-preprocessed_results.R

1. Inputs from pre-existing RDS objects in results/ folder

```
results/xcell_output.rds
results/quantiseq_output.rds
```

2. Running the analysis:

```
bash 03-Immune-deconv-preprocessed_results.sh

```

### Script Changes:

```
#### Script : 03-Immune-deconv-preprocessed_results.R

a. Tow opt parser (filter_cancer_group & filter_for_cancer_subtype) has been added to work on any cancer group and its molecular subtype. In this case, we are filtering for "Medulloblastoma"
b. Plotting for Immune cell distributions for each Cohort (in this case, PBTA, PPTC) for the Cancer molecular subtypes (After filtering, "MB, WNT"; "MB, Group3"; "MB, Group4" and "MB, SHH" were found!)
c. Preliminary summary tables were generated for understanding the pre-processed data (can be ignored!)

```

3. Output: 

```
For (deconv_methods: xcell and quantiseq)

results_custom_deconv_subtyping/{deconv_method}_deconv_summary_table.tsv 
results_custom_deconv_subtyping/{deconv_method}_filtered_output_summary_table.tsv
results_custom_deconv_subtyping/{deconv_method}_immuneCell_distributions_byMolecularSubtype.pdf

```

