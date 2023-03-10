"""
Script to use a GTF and expression matrix with ENSG to liftover gene symbols
"""

import argparse
import gzip
import re
import sys
from statistics import mean


def update_gtf_dict(in_dict, in_text):
    """
    Update ensg - sym ID dictionary on gene entries only
    """
    data = in_text.strip('\n').split('\t')
    if data[tx_type] == 'gene':
        try:
            gene_id = id_find.match(data[tx_desc]).group(1)
            gene_name = name_find.match(data[tx_desc]).group(1)
            in_dict[gene_id] = gene_name
        except Exception as e:
            sys.stderr.write(str(e) + '\n')
            sys.stderr.write("Failed to parse gene ID and name from " + in_text)


def liftover(in_data, out_lift, n, l):
    """
    Update gene symbol if found in dict, leave alone if not
    """
    if l % n == 0:
        sys.stderr.write("Processed " + str(l) + " lines\n")
        sys.stderr.flush()
    data = in_data.rstrip('\n').split('\t')
    cur_id = data[g_idx].split('.')[0]
    if cur_id in gtf_dict:
        data[n_idx] = gtf_dict[cur_id]
    gene_sym = data[n_idx]
    if gene_sym not in sym_ct:
        sym_ct[gene_sym] = 1
    else:
        sym_ct[gene_sym] += 1
    out_data = cur_id + '\t' + gene_sym + '\t' + '\t'.join(data[d_idx:]) + '\n'
    out_lift.write(out_data)
    l += 1
    return l


def collate_dups(in_data, in_dict, ct, out_file, l):
    """
    Checks if entry is in dup dict (really smybol seen > 1 in sym_ct),
    then add to dict ot process, or output  to open file handle as-is if not repeated
    """
    if l % n == 0:
            sys.stderr.write("Processed " + str(l) + " lines\n")
            sys.stderr.flush()
    data = in_data.rstrip('\n').split('\t')
    gene_sym = data[n_idx]
    if sym_ct[gene_sym] > 1:
        if gene_sym not in dup_dict:
            dup_dict[gene_sym] = []
            ct += 1
        # cast as int for later
        dup_dict[gene_sym].append(list(map(float, data[d_idx:])))
    else:
        out_data = gene_sym + "\t" + "\t".join(data[d_idx:]) + '\n'
        out_file.write(out_data)
    l += 1
    return ct, l


def collapse_dups(dup_dict, out_file):
    """
    Get highest mean for repeated gene symbols and output that row
    """
    for gene_sym in dup_dict:
        means = []
        for data in dup_dict[gene_sym]:
            means.append(mean(data))
        top_idx = means.index(max(means))
        out_data = gene_sym + "\t" + "\t".join(list(map(str, dup_dict[gene_sym][top_idx])))
        out_file.write(out_data)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description="Script to use a GTF and expression matrix with ENSG to liftover gene symbols."
    )
    parser.add_argument(
        "-t",
        "--table",
        action="store",
        dest="table",
        help="Input expression matrix. Can be gzipped or plain tsv",
    )
    parser.add_argument(
        "-o", "--output-basename", action="store", dest="out", help="output basename of liftover and collapse if flag given"
    )
    parser.add_argument(
        "-i", "--input-gtf", action="store", dest="gtf", help="GTF file to use as liftover reference. Can be gzipped or plain tsv"
    )
    parser.add_argument(
        "-g", "--gene-id", action="store", dest="gene_id", help="NAME of gene ID (ENSG values) column"
    )
    parser.add_argument(
        "-n", "--gene-name", action="store", dest="gene_name", help="NAME of gene name column, if present to keep if not found"
    )
    parser.add_argument(
        "-s", "--skip", action="store", dest="skip", type=int, help="Number of lines to skip if needed"
    )
    parser.add_argument(
        "-c", "--collapse", action="store_true", dest="collapse", help="If set, will collapse on gene_symbol"
    )

    args = parser.parse_args()

    gtf_dict = {}
    # hardcoded for now, could make opts
    tx_type = 2
    tx_desc = 8
    # set up pattern matcher ahead of time for efficiency
    id_find = re.compile(r"gene_id \"(ENSG\d+)\.\d+\";")
    name_find = re.compile(r".*gene_name \"(\S+)\";")
    sys.stderr.write('Indexing gtf\n')
    with (gzip.open if args.gtf.endswith("gz") else open)(args.gtf, "rt", encoding="utf-8") as gtf:
        for line in gtf:
            if not line.startswith('##'):
                update_gtf_dict(gtf_dict, line)
    # track sym occurrences if collapse flag given for later processing
    sym_ct = {}
    with (gzip.open if args.table.endswith("gz") else open)(args.table, "rt", encoding="utf-8") as table:
        if args.skip:
            sys.stderr.write('Skipping ' + str(args.skip) + ' lines\n')
            for i in range(args.skip):
                skip = next(table)
        head = next(table)
        header = head.rstrip('\n').split('\t')
        # get positions of current gene id and name
        g_idx = header.index(args.gene_id)
        n_idx = header.index(args.gene_name)
        # make assumption that data will start after ID and name
        d_idx = g_idx + 1
        if n_idx > g_idx:
            d_idx = n_idx + 1
        sys.stderr.write('Processing table\n')
        sys.stderr.flush()
        out_lift_fn = args.out + '.liftover.tsv'
        out_lift = open(out_lift_fn, 'w')
        # Some counters to track progress
        n = 1000
        l = 0
        out_head = ('gene_id\tgene_name\t' + '\t'.join(header[d_idx:]))
        out_lift.write(out_head)
        for line in table:
            l = liftover(line, out_lift, n, l)
        # dump out rest of leftover lines
        out_lift.close()
    
    if args.collapse:
        sys.stderr.write("Collapse flag given. Collating dups from liftover\n")
        sys.stderr.flush()
        dup_dict = {}
        dup_ct = 0
        out_collapse_fn = args.out + '.collapsed.tsv'
        out_collapse = open(out_collapse_fn, 'w')
        with (gzip.open if out_lift_fn.endswith("gz") else open)(out_lift_fn, "rt", encoding="utf-8") as lifted:
            # already got header as array, so can skip splitting
            skip = next(lifted)
            out_head = 'gene_name' + '\t' + '\t'.join(header[d_idx:]) + '\n'
            out_collapse.write(out_head)
            l = 0
            for line in lifted:
                dup_ct, l = collate_dups(line, dup_dict, dup_ct, out_collapse, l)
        sys.stderr.write("Processing " + str(dup_ct) + " repeat gene symbols and choosing highest mean expression\n")
        sys.stderr.flush()
        collapse_dups(dup_dict, out_collapse)
        sys.stderr.write("Done!\n")
        sys.stderr.flush()

