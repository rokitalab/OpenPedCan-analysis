"""
Script to use a GTF and expression matrix with ENSG to liftover gene symbols
"""

import argparse
import gzip
import re
import sys
from statistics import mean


def update_gtf_dict(gtf_dict, in_text, tx_type, tx_desc, id_find, name_find):
    """
    Update ensg - sym ID dictionary on gene entries only
    """
    data = in_text.strip('\n').split('\t')
    if data[tx_type] == 'gene':
        try:
            gene_id = id_find.match(data[tx_desc]).group(1)
            gene_name = name_find.match(data[tx_desc]).group(1)
            gtf_dict[gene_id] = gene_name
        except Exception as e:
            print(str(e), file=sys.stderr)
            print("Failed to parse gene ID and name from " + in_text, file=sys.stderr)


def liftover(line_from_infile, gene_sym_ct_dict, g_idx, n_idx, gtf_dict, d_idx):
    """
    Update gene symbol if found in dict, leave alone if not
    """
    data = line_from_infile.rstrip('\n').split('\t')
    cur_id = data[g_idx].split('.')[0]
    if cur_id in gtf_dict:
        data[n_idx] = gtf_dict[cur_id]
    gene_sym = data[n_idx]
    if gene_sym not in gene_sym_ct_dict:
        gene_sym_ct_dict[gene_sym] = 1
    else:
        gene_sym_ct_dict[gene_sym] += 1
    out_data = "{}\t{}\t{}\n".format(cur_id, gene_sym, '\t'.join(data[d_idx:]))
    return out_data


def collate_dups(in_data, dup_gene_sym_dict, ct, gene_sym_ct_dict, d_idx, n_idx):
    """
    Checks if entry is in dup dict (really symbol seen > 1 in gene_sym_ct_dict),
    then add to dict ot process, or output  to open file handle as-is if not repeated
    """
    data = in_data.rstrip('\n').split('\t')
    gene_sym = data[n_idx]
    out_data = None
    if gene_sym_ct_dict[gene_sym] > 1:
        if gene_sym not in dup_gene_sym_dict:
            dup_gene_sym_dict[gene_sym] = []
            ct += 1
        # cast as int for later
        dup_gene_sym_dict[gene_sym].append(list(map(float, data[d_idx:])))
    else:
        out_data = "{}\t{}\n".format(gene_sym, '\t'.join(data[d_idx:]))
    return ct, out_data


def collapse_dups(dup_gene_sym_dict):
    """
    Get highest mean for repeated gene symbols and output that row
    """
    for gene_sym in dup_gene_sym_dict:
        means = [mean(data) for data in dup_gene_sym_dict[gene_sym]]
        top_idx = means.index(max(means))
        out_data = "{}\t{}\n".format(gene_sym, '\t'.join(list(map(str, dup_gene_sym_dict[gene_sym][top_idx]))))
        return out_data


def main():
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
    print('Indexing gtf', file=sys.stderr)
    with (gzip.open if args.gtf.endswith("gz") else open)(args.gtf, "rt", encoding="utf-8") as gtf:
        for line in gtf:
            if not line.startswith('##'):
                update_gtf_dict(gtf_dict, line, tx_type, tx_desc, id_find, name_find)
    # track sym occurrences if collapse flag given for later processing
    gene_sym_ct_dict = {}
    with (gzip.open if args.table.endswith("gz") else open)(args.table, "rt", encoding="utf-8") as table:
        if args.skip:
            print("Skipping {} lines".format(args.skip), file=sys.stderr)
            for i in range(args.skip):
                skip = next(table)
        head = next(table)
        header = head.rstrip('\n').split('\t')
        # get positions of current gene id and name
        g_idx = header.index(args.gene_id)
        n_idx = header.index(args.gene_name)
        # make assumption that data will start after ID and name
        d_idx = max(g_idx, n_idx) + 1
        print('Processing table', file=sys.stderr)
        out_lift_fn = args.out + '.liftover.tsv'
        out_lift = open(out_lift_fn, 'w')
        # Some counters to track progress
        n = 10000
        l = 0
        out_head = "{}\t{}\t{}\n".format('gene_id', 'gene_name', '\t'.join(header[d_idx:]))
        out_lift.write(out_head)
        for line in table:
            if l % n == 0:
                print("Processed {} lines".format(l), file=sys.stderr)
            out_data = liftover(line, gene_sym_ct_dict, g_idx, n_idx, gtf_dict, d_idx)
            out_lift.write(out_data)
            l += 1
        # dump out rest of leftover lines
        out_lift.close()
    
    if args.collapse:
        print("Collapse flag given. Collating dups from liftover", file=sys.stderr)
        dup_gene_sym_dict = {}
        dup_ct = 0
        out_collapse_fn = args.out + '.collapsed.tsv'
        out_collapse = open(out_collapse_fn, 'w')
        with (gzip.open if out_lift_fn.endswith("gz") else open)(out_lift_fn, "rt", encoding="utf-8") as lifted:
            # already got header as array, so can skip splitting
            skip = next(lifted)
            out_head = "{}\t{}\n".format('gene_name', '\t'.join(header[d_idx:]))
            out_collapse.write(out_head)
            l = 0
            for line in lifted:
                dup_ct, out_data = collate_dups(line, dup_gene_sym_dict, dup_ct, gene_sym_ct_dict, d_idx, n_idx)
                if l % n == 0:
                    print("Processed {} lines".format(l), file=sys.stderr)
                if out_data:
                    out_collapse.write(out_data)
                l += 1
        print("Processing {} repeat gene symbols and choosing highest mean expression".format(dup_ct), file=sys.stderr)
        out_data = collapse_dups(dup_gene_sym_dict)
        out_collapse.write(out_data)
        print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
