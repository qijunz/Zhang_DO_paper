#!/usr/bin/env python3

# purpose: extract fasta sequences from ncbi nr database
#          using 1) ncbi.nr.fa
#                2) diamond blastp hits output file

# USAGE: ncbi_nr_seq_extract.py -p nr -q ncbi_nr_search_hsd_query_UniProt.diamond.out -o ncbi_nr_search_out/


import os
import re
import argparse
import subprocess
import pandas as pd

from Bio import SeqIO

def extract_seqs(ncbi_nr, diamond_out, out_dir):
    diamond_out_pd = pd.read_csv(diamond_out, sep="\t", header=None)
    
    query_ids = list(set(diamond_out_pd.loc[:,0].tolist()))
    
    hit_ids = {k:[] for k in query_ids}
    hit_records = {k:[] for k in query_ids}
    
    for query in query_ids:
        hit_ids[query] = diamond_out_pd.loc[diamond_out_pd.loc[:,0] == query,1].tolist()
    
    ncbi_seqs = SeqIO.parse(ncbi_nr, "fasta")
    for ncbi_seq in ncbi_seqs:
        for query in hit_ids:
            if ncbi_seq.id in hit_ids[query]:
                hit_records[query].append(ncbi_seq)
    
    for query in hit_records:
        SeqIO.write(hit_records[query], out_dir+'{}.fa'.format(query), "fasta")


def main(args):
    ncbi_nr = args.ncbi_nr
    diamond_out = args.query_sequence_file
    out_dir = args.output_dir

    extract_seqs(ncbi_nr, diamond_out, out_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-p', '--ncbi_nr',
                        help='NCBI NR fasta file name',
                        type=str,
                        default='')
    parser.add_argument('-q', '--query_sequence_file',
                        help='query sequence file generated from diamond, 1st column is the names of query sequences, 2nd column is the blastp hits seq ids',
                        type=str,
                        default='')             
    parser.add_argument('-o', '--output_dir',
                        help='output directory name for extrated fasta sequences',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)
