## Add formatted sample id column for PedCBio upload

**Module author:** Run Jin ([@runjin326](https://github.com/runjin326)) & Jo Lynne Rokita
Modified by Miguel Brown

Pre-sample name updates include (temporary):
- Adding GTEx group --> `cancer_group` (HISTOLOGY on pedcbio)
- Adding GTEx subgroup --> `harmonized_diagnosis` (CANCER_TYPE_DETAILED on pedcbio)
- Ensuring all tumor samples previously without a `harmonized_diagnosis` have one (`cancer group`, or `cancer group` + `molecular subtype`, in the case of NBL)
- Fix `broad_histology` discrepancy for heme malignancies


Currently, for some of the samples, when multiple DNA or RNA specimens are associated with the same sample, there 
is no column that would distinguish between different aliquots while still tying DNA and RNA together.
This module adds a column called `formatted_sample_id` where the base name is the sample id and additional `tiebreaks` were added when multiple RNA or DNA samples are associated with the same participant.

For PBTA samples, `sample_id` column is used as the basename
- Using `sample_id` column, we can tie all DNA and RNA samples together
- Using `formatted_sample_id` column, we can distinguish amongst multiple DNA or RNA samples 
  - Multiple DNA samples associated with the same sample would use `aliquot_id` as the tie breaker
  - Multiple RNA samples associated with the same sample would use `RNA_library` as the tie breaker 

For TARGET, TCGA, and GTEx samples, `Kids_First_Participant_ID` column is used as the basename
- Using `Kids_First_Participant_ID` column, we can tie all DNA and RNA samples together
- Using `formatted_sample_id` column, we can distinguish amongst multiple DNA or RNA samples 
  - For TARGET, `Kids_First_Participant_ID` + last 7 digits from the `Kids_First_Specimen_ID` is used as formatted sample ID
  - For TCGA, `Kids_First_Participant_ID` + `sample_id` + `aliquot_id` is used as formatted sample ID
  - For GTEx, `Kids_First_Participant_ID` + `aliquot_id` is used as formatted sample ID

For DGD, use `aliquot_id`, but remove test name

Usage:
```sh
# In repo:
Rscript --vanilla pedcbio-sample-name/pedcbio_sample_name_col.R -h path-to-histolgies-file.tsv"

```
or
```sh
# As standanlone in any env with required R packages
Rscript --vanilla pedcbio_sample_name_col.R -i path-to-histolgies-file.tsv -n path-to-cbio-names.csv
```

Input:
- `input/cbio_sample_ids.csv`

Output:
- `results/histologies-formatted-id-added.tsv`

The output files can directly uploaded to S3 buckets for loading into PedCBio or used locally.
