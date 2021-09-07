#!/usr/bin/env python3

# PURPOSE: Reconstruct fastq file from mouse genome removed BAM file
# Why do this:
#   During the pre-process of WGS libraries, in order to remove host (mouse in this case) DNA,
#   the sequenced reads are mapped against host genome (Mus_musculus.GRCm38.98) and get the BAM alignment results.
#   "samtools view -f 4" gives the query sequence itself is unmapped
#   "samtools view -f 8" gives the the mate is unmapped
#   and "samtools view -f 4 -f 8" gives all mouse genome unmapped sequences, and this will be our finial outputs.

# This script take three inputs:
#   1. R1 forward fastq raw reads
#   2. R2 reverse fastq raw reads
#   3. hose genome unmapped SAM file, in which the 4th column are all 0 representing it's unmapped.

# How to use this script (take "Sorcs_WT2" as sample ID):
#   WGS_unmapped_FASTQ_rebuiled_from_mouseDNAremovedBAM.py -s Sorcs_WT2_mouseDNAremoved.sam \
#                                                          -F Sorcs_WT2_R1_trimmomaticTrimmed_paired.fastq \
#                                                          -R Sorcs_WT2_R2_trimmomaticTrimmed_paired.fastq \
#                                                          -Fu Sorcs_WT2_R1_trimmed_paired_mouseDNAremoved.fastq \
#                                                          -Fm Sorcs_WT2_R1_trimmed_paired_mouseDNA.fastq \
#                                                          -Ru Sorcs_WT2_R2_trimmed_paired_mouseDNAremoved.fastq \
#                                                          -Rm Sorcs_WT2_R2_trimmed_paired_mouseDNA.fastq

import argparse
from Bio import SeqIO

def fastq_reconstruction(sam_path, R1_path, R2_path, R1_unmapped_path, R1_mapped_path, R2_unmapped_path, R2_mapped_path):
    # get unmapped seq name set
    unmapped_seq_list = []
    with open(sam_path, "r") as sam_input:
        lines = sam_input.readlines()
        for line in lines:
            if not line.startswith("@"):
                # avoid using set, which require huge RAM, by checking if seq already exist before appending
                if not line.split()[0] in unmapped_seq_list:
                    unmapped_seq_list.append(line.split()[0])

    with open(R1_unmapped_path, "w") as R1_unmapped, open(R1_mapped_path, "w") as R1_mapped, open(R2_unmapped_path, "w") as R2_unmapped, open(R2_mapped_path, "w") as R2_mapped:
        # parsing R1, forward sequences
        with open(R1_path, "r") as R1_input:
            R1_seq_list = list(SeqIO.parse(R1_input, 'fastq'))
            for R1_seq in R1_seq_list:
                if R1_seq.name in unmapped_seq_list:
                    SeqIO.write(R1_seq, R1_unmapped, "fastq")
                else:
                    SeqIO.write(R1_seq, R1_mapped, "fastq")

        # parsing R2, reverse sequences
        with open(R2_path, "r") as R2_input:
            R2_seq_list = list(SeqIO.parse(R2_input, 'fastq'))
            for R2_seq in R2_seq_list:
                if R2_seq.name in unmapped_seq_list:
                    SeqIO.write(R2_seq, R2_unmapped, "fastq")
                else:
                    SeqIO.write(R2_seq, R2_mapped, "fastq")

def main(args):
    sam_path = args.sam_path
    R1_path = args.forward_input
    R2_path = args.reverse_input
    R1_unmapped_path = args.forward_unmapped_output
    R1_mapped_path = args.forward_mapped_output
    R2_unmapped_path = args.reverse_unmapped_output
    R2_mapped_path = args.reverse_mapped_output

    fastq_reconstruction(sam_path, R1_path, R2_path, R1_unmapped_path, R1_mapped_path, R2_unmapped_path, R2_mapped_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-s', '--sam_path',
                        help='SAM file for WGS reads mapping to mouse genome',
                        type=str,
                        default='')
    parser.add_argument('-F', '--forward_input',
                        help='Paired raw forward fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-R', '--reverse_input',
                        help='Paired raw reverse fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-Fu', '--forward_unmapped_output',
                        help='output file for R1 unmapped reads',
                        type=str,
                        default='')
    parser.add_argument('-Fm', '--forward_mapped_output',
                        help='output file for R1 mapped reads',
                        type=str,
                        default='')
    parser.add_argument('-Ru', '--reverse_unmapped_output',
                        help='output file for R2 unmapped reads',
                        type=str,
                        default='')
    parser.add_argument('-Rm', '--reverse_mapped_output',
                        help='output file for R2 mapped reads',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)