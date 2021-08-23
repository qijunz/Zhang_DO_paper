#!/usr/bin/env python3

# PURPOSE: 
# Metagenomic analysis pipeline including:
#   1. de novo assembly
#   2. ORFs prediction
#   3. Taxonomic annotation
#   4. KEGG function annotation
#   5. Meagene abundance quantification
#   6. 

import os
import glob
import math
import datetime
import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO



def main(args):
    R1_path = args.forward_input
    R2_path = args.reverse_input
    out_dir = args.out_dir
    # bin_dir = args.bin_dir
    data_dir = args.data_dir
    sample = args.sample_name

    # setup for individual output folder
    spades_out = out_dir + "spades_out/"
    if not os.path.isdir(spades_out):
        os.mkdir(spades_out)

    prodigal_out = out_dir + "prodigal_out/"
    if not os.path.isdir(prodigal_out):
        os.mkdir(prodigal_out)

    rsem_out = out_dir + "rsem_out/"
    if not os.path.isdir(rsem_out):
        os.mkdir(rsem_out)

    ### de novo assembly, metaSPAdes
    run_spades = "metaspades.py -k 21,33,55,77 \
                                --pe1-1 {} \
                                --pe1-2 {} \
                                -t 4 \
                                -o {}".format(R1_path, R2_path, spades_out)
    subprocess.call(run_spades, stdout=True, shell=True)

    # filter assembled contigs, only keep contigs longer than 500bp
    contig_all = spades_out + "contigs.fasta"

    contig_threshold = 500
    contig_filter = out_dir + "{}.contigs.{}bp.fasta".format(sample, str(contig_threshold))

    seq_len = []
    seq_filter = []

    records = SeqIO.parse(open(contig_all, "r"), "fasta")
    for record in records:
        this_len = len(record.seq)
        seq_len.append(this_len)
        if this_len > contig_threshold:
            seq_filter.append(record)

    SeqIO.write(seq_filter, contig_filter, "fasta")
    print(datetime.datetime.now())
    print('Finish assembly...\n')

    ### ORFs prediction
    run_prodigal = "prodigal -i {} \
                            -a {}/{}.orfs.faa \
                            -d {}/{}.orfs.fna \
                            -p meta -o {}/{}.orfs.out".format(contig_filter, prodigal_out, sample, prodigal_out, sample,prodigal_out, sample)

    subprocess.call(run_prodigal, stdout=True, shell=True)
    
    # filter predicted ORFs, only keep ORFs longer than 100bp
    orf_aa_all = prodigal_out + "{}.orfs.faa".format(sample)
    orf_nt_all = prodigal_out + "{}.orfs.fna".format(sample)

    orf_threshold = 100

    orf_aa_filter = out_dir + "{}.orfs.{}bp.faa".format(sample, str(orf_threshold))
    orf_nt_filter = out_dir + "{}.orfs.{}bp.fna".format(sample, str(orf_threshold))

    # ORF protein translations
    aa_total = 0
    aa_filter = 0
    orf_aa_out = []

    records = SeqIO.parse(open(orf_aa_all, "r"), "fasta")
    for record in records:
        aa_total += 1
        start = record.description.split(" ")[2]
        end = record.description.split(" ")[4]
        record.id = "{}|{}".format(sample, record.id)
        if int(end)-int(start)+1 > orf_threshold:
            aa_filter += 1
            orf_aa_out.append(record)

    SeqIO.write(orf_aa_out, orf_aa_filter, "fasta")

    # ORF nucleotide sequences
    nt_total = 0
    nt_filter = 0
    orf_nt_out = []

    records = SeqIO.parse(open(orf_nt_all, "r"), "fasta")
    for record in records:
        nt_total += 1
        start = record.description.split(" ")[2]
        end = record.description.split(" ")[4]
        record.id = "{}|{}".format(sample, record.id)
        if int(end)-int(start)+1 > orf_threshold:
            nt_filter += 1
            orf_nt_out.append(record)

    SeqIO.write(orf_nt_out, orf_nt_filter, "fasta")

    if aa_total == nt_total and aa_filter == nt_filter:
        print("The number of total predicted ORFs is {}, among which the ORFs with length >100bp is {}.".format(aa_total, aa_filter))
    else:
        print("The number of protein translations with nucleotide sequences is WRONG!")

    print(datetime.datetime.now())
    print('Finish ORFs prediction...\n')

    ### RSEM quantification of metagenes
    # build index
    run_rsem_index = "rsem-prepare-reference -p 4 \
                                             --bowtie2 {} {}{}.rsem.index".format(orf_nt_filter, rsem_out, sample)
    subprocess.call(run_rsem_index, stdout=True, shell=True)

    run_rsem = "rsem-calculate-expression -p 4 \
                                          --bowtie2 \
                                          --paired-end \
                                          --estimate-rspd \
                                          --no-bam-output \
                                          {} \
                                          {} \
                                          {}{}.rsem.index \
                                          {}{}.rsem".format(R1_path, R2_path, rsem_out, sample, rsem_out, sample)
    subprocess.call(run_rsem, stdout=True, shell=True)

    # clear index
    subprocess.call("rm {}/{}.rsem.index*".format(rsem_out, sample), stdout=True, shell=True)
    print(datetime.datetime.now())
    print('Finish RSEM estimation...\n')

    ### KEGG annotation
    # check KOfam data
    ko_profile = data_dir + "/profiles/"
    ko_list = data_dir + "/ko_list"

    if not os.path.isdir(ko_profile):
        print("Please download KOfam profiles!")
    else:
        print("KOfam profiles are ready!")
        
    if not os.path.isfile(ko_list):
        print("Please download KOfam list!")
    else:
        print("KOfam list are ready!")
        

    # make output dir
    kofam_out_file = out_dir + "{}.kofam.out.txt".format(sample)
    
    run_kofam = "exec_annotation -o {} \
                                 --profile {} \
                                 --ko-list {} \
                                 --format mapper-one-line \
                                 --cpu=4 \
                                 {}".format(kofam_out_file, ko_profile, ko_list, orf_aa_filter)

    subprocess.call(run_kofam, stdout=True, shell=True)
    print(datetime.datetime.now())
    print('Finish KO annotation...\n')

    ### Gene abundance sum to KO
    ko_abund_file = out_dir + "{}.gene.tpm.sum2ko.tsv".format(sample)
    rsem_gene_matrix = pd.read_csv(rsem_out + "{}.rsem.genes.results".format(sample), sep="\t")
    gene2ko = pd.read_csv(kofam_out_file, sep="\t", names=["gene", "K1", "K2"])

    rsem_gene_matrix_ko = pd.concat([rsem_gene_matrix[["gene_id", "length", "expected_count", "TPM", "FPKM"]], gene2ko[["K1", "K2"]]], axis=1)
    rsem_gene_matrix_ko.rename(columns={'K1':'KO'}, inplace=True)

    rsem_gene_matrix_ko_sum = rsem_gene_matrix_ko.groupby(["KO"]).sum()[["TPM", "FPKM"]].reset_index()

    ko_list = [x.split("/")[-1].split(".")[0] for x in sorted(glob.glob(ko_profile + '*.hmm'))]
    ko_abund_out = pd.DataFrame(ko_list, columns=['KO'])
    ko_abund_out_merge = pd.merge(ko_abund_out, rsem_gene_matrix_ko_sum, left_on=['KO'], right_on=['KO'], how='left').fillna(0)
    ko_abund_out_merge.to_csv(ko_abund_file, sep="\t", index=False)
    
    print(datetime.datetime.now())
    print('Finish all...\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-F', '--forward_input',
                        help='Paired reads forward fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-R', '--reverse_input',
                        help='Paired reads reverse fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-o', '--out_dir',
                        help='The path of directory for output results',
                        type=str,
                        default='')
    parser.add_argument('-d', '--data_dir',
                        help='The path of directory for (pre)downloaded dataset',
                        type=str,
                        default='')
    parser.add_argument('-b', '--bin_dir',
                        help='The path of directory for bins',
                        type=str,
                        default='')                        
    parser.add_argument('-s', '--sample_name',
                        help='The name of sample used for output file prefix',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)


