#!/usr/bin/env python3


"""
filter-mtp-tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Convert MTP TSV table to JSONL
"""

__author__ = ('Eric Wafula (wafulae@chop.edu)')
__version__ = '1.0'
__date__ = '03 March 2023'


import os 
import sys
import csv
import git
import gzip
import json
import argparse
import numpy as np
import pandas as pd
from collections import OrderedDict

METHLY_VALUES = ["beta", "m"]

def read_parameters():
     p = argparse.ArgumentParser(description=("The 02-mtp-tsv2jsonl.py script converts mtp table from TSV format to JSONL."), formatter_class=argparse.RawTextHelpFormatter)
     p.add_argument('MTP_TSV', type=str, default=None, help="MTP TSV file\n\n")
     p.add_argument('-v', '--version', action='version', version="02-mtp-tsv2jsonl.py script version {} ({})".format(__version__, __date__), help="Print the current 02-mtp-tsv2jsonl.py script version and exit\n\n")
     return p.parse_args()


# establish base dir
root_dir = git.Repo('.', search_parent_directories=True).working_tree_dir

# Set path to module and results directories
module_dir = os.path.join(root_dir, "analyses", "filter-mtp-tables")
mtp_filtered_dir = os.path.join(root_dir, "scratch", "mtp-filtered")


def tsv_to_jsonl(input_tsv_file, output_jsonl_file):
     """Converts MTP table from TSV format to JSONL
     Parameters
     ----------
     input_tsv_file : str
          MTP TSV input file
     output_jsonl_file : str
          MTP JSONL output file name
     Returns
     -------
     None
     """

     tsv_file = gzip.open(input_tsv_file, "rt")
     jsonl_file = gzip.open(output_jsonl_file, "wt")
     reader = csv.DictReader(tsv_file, delimiter="\t")
     headers = reader.fieldnames
     for row in reader:
          row_dict = OrderedDict()
          for header in headers:
               row_dict[header] = row[header]
          json.dump(row_dict, jsonl_file)
          jsonl_file.write("\n")
     tsv_file.close()
     jsonl_file.close()

def main():
     # get input parameter
     args = read_parameters()

     print("========================================")
     print("Converting MTP table from TSV to JSONL")
     print("========================================")
     print("\nConverting MTP table from TSV to JSONL...\n")
     json_file_name = ".".join([os.path.basename(args.MTP_TSV).split(".")[0], "jsonl.gz"])
     tsv_to_jsonl(args.MTP_TSV, "{}".format(json_file_name))

     print("Done converting TSV to JSONL...\n")

if __name__ == "__main__":
     main()