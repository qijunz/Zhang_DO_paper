# purpose: give a list of bins fasta files (take Akkermsnsia as example here)
# and metagenomics reads from one sample, 
# calculate coverage in each bin by provided library
# 
# 
# 

import os
import re
import argparse
import subprocess
import pandas as pd


def run_command(command_string, stdout_path = None):
    # checks if a command ran properly. If it didn't, then print an error message then quit
    print('bins_coverage_quantification.py, run_command: ' + command_string)
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print('bins_coverage_quantification.py: Error, the command:')
        print(command_string)
        print('failed, with exit code ' + str(exit_code))
        exit(1)


def run_bowtie(path2ref, path2f_seq, path2r_seq, num_processors=2, out_dir = None):
    # given reference fasta seq file (un-indexed) and F/R reads file
    # do:
    # 
    # 1. indexing the reference seq
    # e.g., Akk_HQ_bin_DO025_cluster_DBSCAN_round1_2.fasta
    bowtie2_db = path2ref.split('/')[-1].split('.')[0]
    bowtie2_db_path = out_dir + '/' + bowtie2_db
    run_command('bowtie2-build {} {}'.format(path2ref, bowtie2_db_path))

    # 2. align F/R reads to ref seq
    sam_file = out_dir + '/' + bowtie2_db + '.sam'
    run_command('bowtie2 -x {} -1 {} -2 {} --very-sensitive-local --no-unal -p {} -S {}'.format(bowtie2_db_path, path2f_seq, path2r_seq, num_processors, sam_file))

    # 3. convert sam to sorted bam file
    bam_sorted_file = out_dir + '/' + bowtie2_db + '.sorted.bam'
    run_command('samtools view -bS {} | samtools sort -o {}'.format(sam_file, bam_sorted_file))
    # clean up sam file
    run_command('rm {}'.format(sam_file))
    # clean up index files
    [run_command('rm {}/{}'.format(out_dir, f)) for f in os.listdir(out_dir) if re.match('.*bt2$', f)]

    # 4. tabulate the average coverage of each contig
    contig_len_tbl = out_dir + '/' + bowtie2_db + '.length.txt'
    run_command('python fasta_length_table.py -s {} > {}'.format(path2ref, contig_len_tbl))

    # 5. run genomeCoverageBed to get bed file
    bed_file = out_dir + '/' + bowtie2_db + '.bed'
    run_command('genomeCoverageBed -ibam {} -g {} > {}'.format(bam_sorted_file, contig_len_tbl, bed_file))
    # clean up sorted bam file
    run_command('rm {}'.format(bam_sorted_file))

    # 6. building coverage table in this bin
    cov_out = out_dir + '/' + bowtie2_db + '.cov.tab'
    run_command('python contig_coverage_from_bedtools.py -b {} > {}'.format(bed_file, cov_out))

    # return cov_out


def cov_tbl(sample_id, out_dir):
    # after align and get coverage in all bins, integrate a table
    # input: *tab files in out_dir
    # output: sample_id.tsv in same out_dir

    cov_tsv = out_dir + '/' + sample_id + '.cov.tsv'
    with open(cov_tsv, 'w') as out:

        print('bin_name\tcoverage', file = out)
        cov_tbl_list = ['{}/{}'.format(out_dir, f) for f in os.listdir(out_dir) if re.match('.*tab$', f)]

        for cov_tbl in sorted(cov_tbl_list):
            bin_name = cov_tbl.split('/')[-1].split('.')[0]
            cov_pd = pd.read_csv(cov_tbl, sep='\t')
            # take the genome coverage value (coverage on all contigs in bin)
            cov_value = list(cov_pd.coverage)[-1]

            print('{}\t{}'.format(bin_name, cov_value), file=out)


def main(args):

    bin_dir = args.bin_dir
    seq_f = args.forward_reads
    seq_r = args.reverse_reads
    out_dir = args.output_dir

    core = args.processors
    mouse = args.sample_id

    # check if the output dir exists
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # go through all bin fasta files in bin directory
    bin_files = [f for f in os.listdir(bin_dir) if re.match('.*fasta$', f)]

    for bin_file in sorted(bin_files):
        ref = bin_dir + '/' + bin_file
        print(ref)
        run_bowtie(ref, seq_f, seq_r, core, out_dir)

    # get integrated tsv
    cov_tbl(mouse, out_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-b', '--bin_dir',
                        help='Directory path for all bins, each bin is a fasta file',
                        type=str,
                        default='')
    parser.add_argument('-F', '--forward_reads',
                        help='Paired forward fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-R', '--reverse_reads',
                        help='Paired reverse fastq file (*.fq/*.fastq)',
                        type=str,
                        default='')
    parser.add_argument('-s', '--sample_id',
                        help='ID for fastq seq sample',
                        type=str,
                        default='')
    parser.add_argument('-p', '--processors',
                        help='Number of CPU',
                        type=int,
                        default=1)
    parser.add_argument('-o', '--output_dir',
                        help='output directory contains output files including tsv file with bin name in 1st column and coverage in 2nd column',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)

