---
output:
  html_document: default
  pdf_document: default
---
# Data Formats in Data Download

The release notes for each release are provided in the `release-notes.md` file that accompanies the data files.
A table with brief descriptions for each data file is provided in the `data-files-description.md` file included in the download.

## Processed Data Files

Processed data files are all files derived from samples (e.g., tumors, cell lines) that are processed upstream of this repository and are not the product of any analysis code in the `AlexsLemonade/OpenPBTA-analysis` or `PediatricOpenTargets/OpenPedCan-analysis` repository.

### Consensus Somatic Variant Data
Somatic calls that are retained if they are supported by at least 2 callers OR marked as `HotSpotAllele` because they overlap SNV/INDELs considered as [Cancer Hotspots](https://www.cancerhotspots.org/#/download) OR are TERT promoter SNVs. Please find additional information [here](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md)

* `snv-consensus-plus-hotspots.maf.tsv.gz`

### Somatic Copy Number Variant (CNV) Data

Somatic Copy Number Variant (CNV) data are provided in a modified [SEG format](https://software.broadinstitute.org/software/igv/SEG) for each of the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#somatic-copy-number-variant-calling) and denoted with the `cnv` prefix.
**Somatic copy number data is only generated for whole genome sequencing (WGS) samples.**

* `cnv-cnvkit.seg.gz` is the the CNVkit SEG file. This file contains an additional column `copy.num` to denote copy number of each segment, derived from the CNS file output of the algorithm [described here](https://cnvkit.readthedocs.io/en/stable/fileformats.html).
* `cnv-controlfreec.tsv.gz` is the ControlFreeC TSV file. It is a merge of `*_CNVs` files produced from the algorithm, and columns are [described here](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT).
  
#### A Note on Ploidy

The _copy number_ annotated in the CNVkit SEG file is annotated with respect to ploidy 2, however, the _status_ annotated in the ControlFreeC TSV file is annotated with respect to inferred ploidy from the algorithm, which is recorded in the `histologies.tsv` file. 

### Gene Expression Data

Gene expression estimates from the [applied software packages](https://alexslemonade.github.io/OpenPBTA-manuscript/#gene-expression-abundance-estimation) are provided as a feature (e.g., gene or transcript) by sample matrix.
Gene expression are available in multiple forms in the following files:

* `gene-counts-rsem-expected_count.rds`
* `gene-expression-rsem-tpm.rds`

See [the data description file](data-description.md) for more information about the individual gene expression files.


If your analysis requires de-duplicated gene symbols as row names, please use the collapsed matrices provided as part of the data download ([see below](#collapsed-expression-matrices)).

### RNA splice events Data

RNA splice events rmats file generated from the [established pipeline](https://github.com/d3b-center/pbta-splicing/blob/main/workflows/rmats_wf.cwl) is provided as:

  * `splice-events-rmats.tsv.gz`

### Derived Fusion Files

The filtered and prioritized fusion and downstream files are a product of the [`analyses/fusion_filtering`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion_filtering) analysis module. 

  * `fusion-putative-oncogenic.tsv` contains the filtered and prioritized fusions. 

Binary matrices for the presence of tumor-specific fusions across all RNA biospecimens are the product of [`fusion-summary`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion-summary).

  * `fusion_summary_embryonal_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent embryonal tumor fusions of interest per individual RNA-seq specimen.
  * `fusion_summary_ependymoma_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent ependymal tumor fusions of interest per individual RNA-seq specimen.
  * `fusion_summary_ewings_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent Ewing's sarcoma tumor fusions of interest per individual RNA-seq specimen.
  * `fusion_summary_lgat_foi.tsv` contains a binary matrix that denotes the presence or absence of a recurrent LGAT tumor fusions of interest per individual RNA-seq specimen.


### Structural Variant Data

Structural Variants data produced by the [`MANTA` package](https://github.com/Illumina/manta) is 
provided as:

* `sv-manta.tsv.gz`

### Proteomic data

Whole cell proteomic and phosphorylation proteomic data from `project HOPE` and `CPTAC` are provided as: 

| File name | Data type | Data source | Data description |
|---|---|---|---|
|`cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz`| Processed data | [CPTAC pediatric brain tumor phospho-proteomics expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed phospho-protein expression, log2 abundance
|`cptac-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | [CPTAC pediatric brain tumor protein expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed whole cell protein expression, total abundance
|`cptac-protein-imputed-prot-expression-log2-ratio.tsv.gz`| Processed data | [CPTAC pediatric brain tumor protein expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed whole cell protein expression, log2 abundance
|`gbm-protein-imputed-phospho-expression-abundance.tsv.gz`| Processed data | [CPTAC adult GBM brain tumor phospho-proteomics expression](https://www.sciencedirect.com/science/article/pii/S1535610821000507) | Imputed phospho-protein expression, total abundance
|`gbm-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | [CPTAC adult GBM brain tumor protein expression](https://www.sciencedirect.com/science/article/pii/S1535610821000507) | Imputed whole cell expression, total abundance
|`hope-protein-imputed-phospho-expression-abundance.tsv.gz`| Processed data | Adult and Young Adolescent (AYA) brain tumor phospho-proteomics expression (Project HOPE) | Imputed phospho-protein expression, total abundance
|`hope-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | Adult and Young Adolescent (AYA) brain tumor protein expression (Project HOPE) | Imputed whole cell protein expression, total abundance

### Harmonized Clinical Data

[Harmonized clinical data](https://alexslemonade.github.io/OpenPBTA-manuscript/#clinical-data-harmonization) are released as tab separated values in the following files:

* `histologies.tsv`
* `histologies-base.tsv`

### Independent Sample Lists

[Independent sample list](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/independent-samples) are released as tab separated values in the following files. 
`wgswxspanel` indicates it includes all experimental strategies for DNA sequencing, `rnaseqpanel` indicates it includes all experimental strategies for RNA sequencing.
`eachcohort` indicates the selection is cohort-based. 
Additionally, `.prefer.wxs` (or `.prefer.wgs`) indicates WXS (or WGS) samples were preferentially select when both WGS and WXS are available for a particular participant. 

* `independent-specimens.methyl.primary-plus.tsv`
* `independent-specimens.methyl.primary.tsv`
* `independent-specimens.methyl.relapse.tsv`
* `independent-specimens.rnaseq.primary-plus-pre-release.tsv`
* `independent-specimens.rnaseq.primary-pre-release.tsv`
* `independent-specimens.rnaseq.primary.eachcohort.tsv`
* `independent-specimens.rnaseq.primary.tsv`
* `independent-specimens.rnaseq.relapse-pre-release.tsv`
* `independent-specimens.rnaseq.relapse.eachcohort.tsv`
* `independent-specimens.rnaseq.relapse.tsv`
* `independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv`
* `independent-specimens.rnaseqpanel.primary-plus.pre-release.tsv`
* `independent-specimens.rnaseqpanel.primary-plus.tsv`
* `independent-specimens.rnaseqpanel.primary.eachcohort.tsv`
* `independent-specimens.rnaseqpanel.primary.pre-release.tsv`
* `independent-specimens.rnaseqpanel.primary.tsv`
* `independent-specimens.rnaseqpanel.relapse.eachcohort.tsv`
* `independent-specimens.rnaseqpanel.relapse.pre-release.tsv`
* `independent-specimens.rnaseqpanel.relapse.tsv`
* `independent-specimens.wgs.primary-plus.eachcohort.tsv`
* `independent-specimens.wgs.primary-plus.tsv`
* `independent-specimens.wgs.primary.eachcohort.tsv`
* `independent-specimens.wgs.primary.tsv`
* `independent-specimens.wgs.relapse.eachcohort.tsv`
* `independent-specimens.wgs.relapse.tsv`
* `independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.primary.eachcohort.tsv`
* `independent-specimens.wgswxspanel.primary.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.primary.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.primary.tsv`
* `independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.relapse.eachcohort.tsv`
* `independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv`
* `independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv`
* `independent-specimens.wgswxspanel.relapse.tsv`

## Analysis Files

Analysis files are created by a script in `analyses/*`. 
They can be viewed as _derivatives_ of Processed data files.

### Collapsed Expression Matrices

Collapsed expression matrices are products of the [`analyses/collapse-rnaseq`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq) analysis module.
In cases where more than one Ensembl gene identifier maps to the same gene symbol, the instance of the gene symbol with the maximum mean FPKM in the RSEM FPKM file is retained to produce the following files:

* `gene-counts-rsem-expected_count-collapsed.rds`  
* `gene-expression-rsem-tpm-collapsed.rds`
* `rna-isoform-expression-rsem-tpm.rds`

Additionally, available TCGA and GTEx gene expression files with same format are included:
* `tcga-gene-counts-rsem-expected_count.rds`
* `tcga-gene-expression-rsem-tpm.rds`
* `gtex_gene-counts-rsem-expected_count-collapsed.rds`
* `gtex_gene-expression-rsem-tpm-collapsed.rds`

### Derived Copy Number Files

#### Consensus Copy Number File

Copy number consensus calls from the copy number and structural variant callers are a product of the [`analyses/copy_number_consensus_call`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call) analysis module. 

* `cnv-consensus.seg.gz` contains consensus segments and segment means (log R ratios) from two or more callers, as described in the [analysis README](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call/README.md) - contains only WGS samples. 

* `cnvkit_with_status.tsv` and `consensus_seg_with_status.tsv` contain CNVkit calls for WXS or CNV consensus calls for WGS with gain/loss status, respectively

##### Focal Copy Number Files

Focal copy number files map the consensus calls (genomic segments) in WGS samples to genes for downstream analysis and are a product of the [`analysis/focal-cn-file-preparation`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation).
Note: these files contain biospecimens and genes with copy number changes. 

  - `consensus_seg_annotated_cn_autosomes.tsv.gz` contains focal gene copy number alterations for all autosomes.
  - `consensus_seg_annotated_cn_x_and_y.tsv.gz` contains focal gene copy number alterations for the sex chromosomes.
  
  Focal copy number files in WXS samples only uses results from CNVkit and no consensus calling is required. [`analysis/focal-cn-file-preparation`](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation).
Note: these files contain biospecimens and genes with copy number changes. 

  - `cnvkit_annotated_cn_wxs_autosomes.tsv.gz` contains focal gene copy number alterations for all autosomes.
  - `cnvkit_annotated_cn_wxs_x_and_y.tsv.gz` contains focal gene copy number alterations for the sex chromosomes.
  
  Additionally, autosomes file and x_and_y file for either WGS or WXS to generate the two combined files as followed:
  - `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz`
  - `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_x_and_y.tsv.gz`

  And these two files are further merged to generate:
  - `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz`
