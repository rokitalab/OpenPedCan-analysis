#!/bin/bash
python3 liftover_collapse_rnaseq.py $@ && ls *.tsv | xargs -IFN -P 2 gzip FN