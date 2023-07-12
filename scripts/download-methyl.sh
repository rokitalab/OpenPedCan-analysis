#!/bin/bash

set -e
set -o pipefail

# download methylation files from the OpenPedCan data v12 data release s3 bucket

URL="https://d3b-openaccess-us-east-1-prd-pbta.s3.amazonaws.com/open-targets"
RELEASE="v13"

# check if the release folder exists, if not, create a release folder
[ ! -h "../data/$RELEASE/" ] && mkdir ../data/$RELEASE/

FILES=(`tr -s ' ' < ../data/$RELEASE/methyl-md5sum.txt | cut -d ' ' -f 2` release-notes.md)
for file in "${FILES[@]}"
do
  if [ ! -e "../data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl $URL/$RELEASE/$file -o ../data/$RELEASE/$file
  fi
done

#check md5sum
cd ../data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c methyl-md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done
