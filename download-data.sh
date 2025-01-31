#!/bin/bash

set -e
set -o pipefail

# Use the OpenPedCan bucket as the default.
URL=${OPENPEDCAN_URL:-https://s3.amazonaws.com/bti-openaccess-us-east-1-prd-opc}
RELEASE=${OPENPEDCAN_RELEASE:-v15}
PREVIOUS=${OPENPEDCAN_RELEASE:-v14}

# Remove old symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

# Download the items in FILES
for file in "${FILES[@]}"
do
  if [ ! -e "data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl $URL/$RELEASE/$file -o data/$RELEASE/$file
  fi
done

# Download reference and gencode files from public ftp if do not already exist
GENCODE27="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz"
cd data
if [ ! -e ${GENCODE27##*/} ]
then
  echo "Downloading ${GENCODE27##*/}"
  curl -O $GENCODE27
fi

GENCODE38="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
if [ ! -e ${GENCODE38##*/} ]
then
  echo "Downloading ${GENCODE38##*/}"
  curl -O $GENCODE38
fi


GENCODE39="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz"
if [ ! -e ${GENCODE39##*/} ]
then
  echo "Downloading ${GENCODE39##*/}"
  curl -O $GENCODE39
fi


# if in CI, then we want to generate the reference FASTA from the BSgenome.Hsapiens.UCSC.hg38 R package
# because it is considerably faster to do so

if [ "$RELEASE" == "testing" ]; then
  Rscript -e "Biostrings::writeXStringSet(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, c('chr21', 'chr22', 'chrX', 'chrY')), 'GRCh38.primary_assembly.genome.fa.gz', format = 'fasta', compress = 'gzip')"
else
  REFERENCE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz"
  if [ ! -e ${REFERENCE##*/} ]
  then
    echo "Downloading ${REFERENCE##*/}"
    curl -O $REFERENCE
  fi
fi
cd -

# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# make data directory unwritable in CI
if [ "$RELEASE" == "testing" ]; then
  chmod u-w data
fi
