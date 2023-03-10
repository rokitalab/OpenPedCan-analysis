#!/bin/bash
python3 liftover_collapse_rnaseq.py $@ && gzip *.tsv