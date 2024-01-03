## Data file descriptions

This document contains information about all data files associated with this project. Each file will have the following association information:

+ **File type** will be one of:
	+ *Reference file*: Obtained from an external source/database. When known, the obtained data and a link to the external source is included.
	+ *Modified reference file*: Obtained from an external source/database but modified for OpenPBTA use. 
	+ *Processed data file*: Data that are processed upstream of the analysis project, e.g., the output of a somatic single nucleotide variant method. Links to the relevant D3B Center or Kids First workflow (and version where applicable) are included in **Origin**.
	+ *Analysis file*: Any file created by a script in `analyses/*`. 
+ **Origin**
	+ For _Processed data files_, a link the relevant D3B Center or Kids First workflow (and version where applicable).
	+ When applicable, a link to the specific *script* that produced (or modified, for *Modified reference file* types) the data.
+ **File description**
	+ A *brief* one sentence description of what the file contains (e.g., bed files contain coordinates for features XYZ).

### current release (v13)
| **File name** |  **File Type** | **Origin** | **File Description** |
|---------------|----------------|------------------------|-----------------------|
|`histologies-base.tsv` | Data file | Cohort-specific data files and databases | Clinical and sequencing metadata for each biospecimen
|`histologies.tsv` | Modified data file | [`molecular-subtyping-integrate`](https://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/molecular-subtyping-integrate) | `histologies-base.tsv` plus `molecular_subtype`, `cancer_group`, `integrated_diagnosis`, and `harmonized_diagnosis`
|`intersect_cds_lancet_strelka_mutect_WGS.bed` | Analysis file | [`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/) | Intersection of `gencode.v39.primary_assembly.annotation.gtf.gz` CDS with Lancet, Strelka2, Mutect2 regions
|`intersect_strelka_mutect_WGS.bed` | Analysis file | [`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/) | Intersection of `gencode.v39.primary_assembly.annotation.gtf.gz` CDS with Strelka2 and Mutect2 regions called
|`efo-mondo-map.tsv` | Reference mapping file | [Manual collation](https://github.com/d3b-center/ticket-tracker/issues/88) | Mapping of EFO and MONDO codes to cancer groups
|`efo-mondo-map-prefill.tsv` | Modified reference mapping file | Analysis file generated in `molecular-subtyping-integrate` | Mapping of EFO and MONDO codes to cancer groups
|`ensg-hugo-pmtl-mapping.tsv` | Reference mapping file | [Manual curation of PMTLv1.1 by FNL; RNA-Seq pipeline GTF mapping](https://github.com/d3b-center/ticket-tracker/issues/125) | File which maps Hugo Symbols to ENSEMBL gene IDs an each ENSG to the RMTL curated by FNL
|`*.bed` | Reference file | [Manual collation](https://github.com/d3b-center/ticket-tracker/issues/258) | Bed files used for variant calling and are used for tmb calculation
|`uberon-map-gtex-group.tsv` | Reference mapping file | [Manual collation](https://github.com/d3b-center/ticket-tracker/issues/85) | Mapping of UBERON codes to tissue types in GTEx broad groups
|`uberon-map-gtex-subgroup.tsv` | Reference mapping file | [Manual collation](https://github.com/d3b-center/ticket-tracker/issues/85) | Mapping of UBERON codes to tissue types in GTEx subgroups
|`methyl-beta-values.rds` | Processed data file | [methylation beta values](https://github.com/d3b-center/ticket-tracker/issues/269)| Methylation beta values
|`methyl-m-values.rds` | Processed data file | [methylation m values](https://github.com/d3b-center/ticket-tracker/issues/269)| Methylation m values
|`rna-isoform-expression-rsem-tpm.rds` | Processed data file | [RNA isoform TPM files](https://github.com/d3b-center/ticket-tracker/issues/341)| RNA isoform TPM files
|`fusion-dgd.tsv` | Processed data file | [DGD merged fusion results](https://github.com/d3b-center/ticket-tracker/issues/249)| DGD merged fusion results
|`fusion-arriba.tsv.gz` | Processed data file | [Gene fusion detection](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#gene-fusion-detection); [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_RNAseq_workflow.cwl) | Fusion - [Arriba TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/arriba-tsv-header.md), annotated with FusionAnnotator
|`fusion-starfusion.tsv.gz` | Processed data file | [Gene fusion detection](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#gene-fusion-detection); [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_RNAseq_workflow.cwl) | Fusion - [STARFusion TSV](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/format/starfusion-tsv-header.md)
|`fusion-annoFuse.tsv.gz` | Processed data file | AnnoFuse QC filtered fusion file; [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_annoFuse_wf.cwl) | Filter out normal and non-expressed fusions |
|`fusion_summary_embryonal_foi.tsv` | Analysis file | [`fusion-summary`](ttps://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/fusion-summary) | Summary file for presence of embryonal tumor fusions of interest |
|`fusion_summary_ependymoma_foi.tsv` | Analysis file | [`fusion-summary`](ttps://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/fusion-summary) | Summary file for presence of ependymal tumor fusions of interest |
|`fusion_summary_ewings_foi.tsv` | Analysis file | [`fusion-summary`](ttps://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/fusion-summary) | Summary file for presence of Ewing's sarcoma fusions of interest |
|`fusion_summary_lgg_hgg_foi.tsv` | Analysis file | [`fusion-summary`](ttps://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/fusion-summary) | Summary file for presence of LGG and HGG fusions of interest |
|`fusion-putative-oncogenic.tsv` | Analysis file | [`fusion_filtering`](https://github.com/d3b-center/OpenPedCan-analysis/tree/master/analyses/fusion_filtering) | Filtered and prioritized fusions
|`gene-counts-rsem-expected_count-collapsed.rds` | Analysis file | PBTA+GMKF+TARGET [`collapse-rnaseq`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq) | Gene expression - RSEM expected_count for each samples collapsed to gene symbol (gene-level)
|`gene-expression-rsem-tpm-collapsed.rds` | Analysis file | PBTA+GMKF+TARGET [`collapse-rnaseq`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/collapse-rnaseq) | Gene expression - RSEM TPM for each samples collapsed to gene symbol (gene-level)
|`tcga-gene-expression-rsem-tpm-collapsed.rds` | Modified reference file | TCGA samples lifted from GENCODE v27 to v39 | Gene expression - RSEM TPM for each samples collapsed to gene symbol (gene-level)
|`gtex-gene-expression-rsem-tpm-collapsed.rds` | Modified reference file | GTEX v8 release lifted to GENCODE v39 | Gene expression - RSEM TPM for each samples collapsed to gene symbol (gene-level)
|`gtex_gene-counts-rsem-expected_count-collapsed.rds` | Modified reference file | GTEX v8 release lifted to GENCODE v39 | Gene expression - RSEM counts for each samples collapsed to gene symbol (gene-level)
|`WGS.hg38.lancet.300bp_padded.bed` | Reference Target/Baits File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | WGS.hg38.lancet.unpadded.bed file with each region padded by 300 bp
|`WGS.hg38.lancet.unpadded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) |  hg38 WGS regions created using UTR, exome, and start/stop codon features of the GENCODE 31 reference, augmented with PASS variant calls from Strelka2 and Mutect2
|`WGS.hg38.mutect2.vardict.unpadded.bed` | Reference Regions File  | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) |  hg38 BROAD Institute interval calling list (restricted to Chr1-22,X,Y,M and non-N regions) used for Mutect2 and VarDict variant callers
|`WGS.hg38.strelka2.unpadded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | hg38 BROAD Institute interval calling list (restricted to Chr1-22,X,Y,M) used for Strelka2 variant caller
|`WGS.hg38.vardict.100bp_padded.bed` | Reference Regions File | [SNV and INDEL calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#snv-and-indel-calling) | `WGS.hg38.mutect2.vardict.unpadded.bed` with each region padded by 100 bp used for VarDict variant caller
|`snv-consensus-plus-hotspots.maf.tsv.gz` | Analysis file | [Kids First somatic workflow consensus calls](https://github.com/kids-first/kf-somatic-workflow/blob/master/docs/kfdrc-consensus-calling.md) | Consensus (2 of 4) maf +1/4 hotspots 
|`snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz` | Analysis file | [Kids First Tumor Only workflow](https://github.com/kids-first/kf-tumor-workflow/) | Mutect2 tumor only with additional filters to remove t_alt_count <5
|`cnv-cnvkit.seg.gz` | Processed data file | [Copy number variant calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-copy-number-variant-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc_combined_somatic_wgs_cnv_wf.cwl) | Somatic Copy Number Variant - CNVkit [SEG file](https://cnvkit.readthedocs.io/en/stable/fileformats.html#seg)
|`cnv-consensus.seg.gz` | Analysis file | [`copy_number_consensus_call`]](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call) | Somatic Copy Number Variant - WGS samples only
|`cnvkit_with_status.tsv` <br> `consensus_seg_with_status.tsv` | Analysis files | [`copy_number_consensus_call`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call) | CNVkit calls for WXS or CNV consensus calls for WGS with gain/loss status
|`cnv-consensus-gistic.gz` | Analysis file | [`run-gistic`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/run-gistic) | GISTIC results - WGS samples only
|`cnv-controlfreec.tsv.gz` | Processed data file | [Copy number variant calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-copy-number-variant-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc_combined_somatic_wgs_cnv_wf.cwl) | Somatic Copy Number Variant - TSV file that is a merge of [ControlFreeC `*_CNVs` files](http://boevalab.inf.ethz.ch/FREEC/tutorial.html#OUTPUT)
|`cnv-gatk.seg.gz` | Processed data file | Copy number variant calling | Somatic Copy Number Variant - TSV SEG file produced by GATK CNV | Analysis file | 
|`consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz` | Analysis file | [`focal-cn-file-preparation`](hhttps://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation) | TSV file containing genes with copy number changes per biospecimen; all chromosomes
|`consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_x_and_y.tsv.gz` | Analysis file | [`focal-cn-file-preparation`](hhttps://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation) | TSV file containing genes with copy number changes per biospecimen; sex chromosomes only
|`consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz` | Analysis file | [`focal-cn-file-preparation`](hhttps://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation) | TSV file containing genes with copy number changes per biospecimen; autosomes chromosomes only
|`snv-mutation-tmb-all.tsv` | Processed data file | [`tmb-calculation`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation) | TSV file with sample names and their tumor mutation burden counting all variants
|`snv-mutation-tmb-coding.tsv` | Processed data file | [`tmb-calculation`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation) | TSV file with sample names and their tumor mutation burden counting all variants in coding region only
|`sv-manta.tsv.gz`| Processed data file | [Structural variant calling](https://github.com/AlexsLemonade/OpenPBTA-manuscript/blob/master/content/03.methods.md#somatic-structural-variant-calling); [Workflow](https://github.com/d3b-center/OpenPBTA-workflows/blob/master/cwl/kfdrc_strelka2_mutect2_manta_workflow.cwl) | Somatic Structural Variant - Manta output, annotated with AnnotSV (WGS samples only)
|`splice-events-rmats.tsv.gz`| Processed data | [Kids First splice variant workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/workflow/kfdrc_RNAseq_workflow.cwl); [Workflow](https://github.com/kids-first/kf-rnaseq-workflow/blob/master/tools/rmats_both_bam.cwl) | rMATs single sample workflow
|`cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz`| Processed data | [CPTAC pediatric brain tumor phospho-proteomics expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed phospho-protein expression, log2 abundance
|`cptac-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | [CPTAC pediatric brain tumor protein expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed whole cell protein expression, total abundance
|`cptac-protein-imputed-prot-expression-log2-ratio.tsv.gz`| Processed data | [CPTAC pediatric brain tumor protein expression](https://pubmed.ncbi.nlm.nih.gov/33242424/) | Imputed whole cell protein expression, log2 abundance
|`gbm-protein-imputed-phospho-expression-abundance.tsv.gz`| Processed data | [CPTAC adult GBM brain tumor phospho-proteomics expression](https://www.sciencedirect.com/science/article/pii/S1535610821000507) | Imputed phospho-protein expression, total abundance
|`gbm-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | [CPTAC adult GBM brain tumor protein expression](https://www.sciencedirect.com/science/article/pii/S1535610821000507) | Imputed whole cell expression, total abundance
|`hope-protein-imputed-phospho-expression-abundance.tsv.gz`| Processed data | Adult and Young Adolescent (AYA) brain tumor phospho-proteomics expression (Project HOPE) | Imputed phospho-protein expression, total abundance
|`hope-protein-imputed-prot-expression-abundance.tsv.gz`| Processed data | Adult and Young Adolescent (AYA) brain tumor protein expression (Project HOPE) | Imputed whole cell protein expression, total abundance
|`rna-dna-qc-stats.tsv`| Reference QC file | Quality control metrics for WGS, WXS, DNA panel, and RNA-Seq samples | Used to filter samples for data release
`independent-specimens.methyl.primary.tsv` <br>
`independent-specimens.methyl.relapse.tsv` <br>
`independent-specimens.rnaseq.primary.eachcohort.tsv` <br>
`independent-specimens.rnaseq.primary.tsv` <br>
`independent-specimens.rnaseq.relapse-pre-release.tsv` <br>
`independent-specimens.rnaseq.relapse.eachcohort.tsv` <br>
`independent-specimens.rnaseq.relapse.tsv` <br>
`independent-specimens.rnaseq.primary-plus-pre-release.tsv` <br>
`independent-specimens.rnaseqpanel.primary-plus.pre-release.tsv` <br>
`independent-specimens.rnaseqpanel.primary-plus.tsv` <br>
`independent-specimens.rnaseqpanel.primary.eachcohort.tsv` <br>
`independent-specimens.rnaseqpanel.primary.tsv` <br>
`independent-specimens.rnaseqpanel.relapse.eachcohort.tsv` <br>
`independent-specimens.rnaseqpanel.relapse.tsv` <br>
`independent-specimens.wgs.primary-plus.eachcohort.tsv` <br>
`independent-specimens.wgs.primary-plus.tsv` <br>
`independent-specimens.wgs.primary.eachcohort.tsv` <br>
`independent-specimens.wgs.primary.tsv` <br>
`independent-specimens.wgs.relapse.eachcohort.tsv` <br>
`independent-specimens.wgs.relapse.tsv` <br>
`independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.primary.eachcohort.tsv` <br>
`independent-specimens.wgswxspanel.primary.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.primary.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.primary.tsv` <br>
`independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.relapse.eachcohort.tsv` <br>
`independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv` <br>
`independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv` <br>
`independent-specimens.wgswxspanel.relapse.tsv` <br>
`independent-specimens.methyl.primary-plus.eachcohort.tsv` <br>
`independent-specimens.methyl.primary.eachcohort.tsv` <br>
`independent-specimens.methyl.relapse.eachcohort.tsv`| Analysis files | [`independent-samples`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/independent-samples) | Independent (non-redundant) sample list of DNA, RNA, or methylation samples of all sequencing methods, from primary, primary-plus, or relapse tumors within each or across all cohorts
`independent-specimens.rnaseq.primary-plus-pre-release.tsv` <br>
`independent-specimens.rnaseq.primary-pre-release.tsv` <br>
`independent-specimens.rnaseq.primary-pre-release.tsv` <br>
`independent-specimens.rnaseq.relapse-pre-release.tsv` | Analysis files | [`independent-samples`](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/independent-samples) | Independent (non-redundant) sample list of RNA samples of all sequencing methods, from primary, primary-plus, or relapse tumors across all cohorts for the purposes of running `fusion_filtering` pre-release
