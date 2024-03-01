#!/bin/bash

# R. Corbett (adapted from J. Taroni for ALSF CCDL)
# December 2022

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"


# In CI we'll run an abbreviated version of the de novo signatures extraction
ABBREVIATED_MUTSIGS=${OPC_QUICK_MUTSIGS:-0}

# Run only consensus testing file in CI, since tumor only snv is large
IS_CI=${OPENPBTA_TESTING:-0}

if [[ "$IS_CI" -eq "1" ]]

then

  echo "Run the SBS mutational signatures analysis using existing signatures on consensus SNV"
  Rscript -e "rmarkdown::render('01-known_signatures.Rmd', params = list(snv_file = \"snv-consensus-plus-hotspots.maf.tsv.gz\", output_Folder = \"ConsensusSNV\"), clean = TRUE)"
  mv 01-known_signatures.nb.html 01-ConsensusSNV_known_signatures.nb.html

  echo "Run the mutational signatures analysis using COSMIC DBS signatures (v3.3) on consensus SNV"
  Rscript -e "rmarkdown::render('02-cosmic_dbs_signatures.Rmd', params = list(snv_file = \"snv-consensus-plus-hotspots.maf.tsv.gz\", output_Folder = \"ConsensusSNV\"), clean = TRUE)"
  mv 02-cosmic_dbs_signatures.nb.html 02-ConsensusSNV_cosmic_dbs_signatures.nb.html

  echo "Run analysis of adult CNS mutational signatures on consensus SNV"
  Rscript --vanilla 03-fit_cns_signatures.R \
    --snv_file snv-consensus-plus-hotspots.maf.tsv.gz \
    --output_Folder ConsensusSNV 

else 

  # Run the SBS mutational signatures analysis using existing signatures on consensus SNV
  Rscript -e "rmarkdown::render('01-known_signatures.Rmd', params = list(snv_file = \"snv-consensus-plus-hotspots.maf.tsv.gz\", output_Folder = \"ConsensusSNV\"), clean = TRUE)"
  mv 01-known_signatures.nb.html 01-ConsensusSNV_known_signatures.nb.html
  
  # Run the SBS mutational signatures analysis using existing signatures on tumor only SNV
  Rscript -e "rmarkdown::render('01-known_signatures.Rmd', params = list(snv_file = \"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz\", output_Folder = \"TumorOnlySNV\"), clean = TRUE)"
  mv 01-known_signatures.nb.html 01-TumorOnly_known_signatures.nb.html

  # Run the mutational signatures analysis using COSMIC DBS signatures (v3.3)
  Rscript -e "rmarkdown::render('02-cosmic_dbs_signatures.Rmd', params = list(snv_file = \"snv-consensus-plus-hotspots.maf.tsv.gz\", output_Folder = \"ConsensusSNV\"), clean = TRUE)"
  mv 02-cosmic_dbs_signatures.nb.html 02-ConsensusSNV_cosmic_dbs_signatures.nb.html

  Rscript -e "rmarkdown::render('02-cosmic_dbs_signatures.Rmd', params = list(snv_file = \"snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz\", output_Folder = \"TumorOnlySNV\"), clean = TRUE)"
  mv 02-cosmic_dbs_signatures.nb.html 02-TumorOnly_cosmic_dbs_signatures.nb.html
  
  # Run analysis of adult CNS mutational signatures
  Rscript --vanilla 03-fit_cns_signatures.R \
    --snv_file snv-consensus-plus-hotspots.maf.tsv.gz \
    --output_Folder ConsensusSNV 

  Rscript --vanilla 03-fit_cns_signatures.R \
    --snv_file snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz \
    --output_Folder TumorOnlySNV 

  Rscript -e "rmarkdown::render('04-explore_hypermutators.Rmd', params = list(output_Folder = \"ConsensusSNV\"), clean = TRUE)"

  ## Tumor only result did not have sample passing the filter. Therefore, 04 script is not running for tumor only
  #Rscript -e "rmarkdown::render('04-explore_hypermutators.Rmd', params = list(output_Folder = \"TumorOnlySNV\"), clean = TRUE)"

fi
