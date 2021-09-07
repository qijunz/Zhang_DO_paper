# Program to print a table of coverage from bed file results
# USAGE: contig_coverage_from_bedtools.py -b <input bedfile> [> out]
# prints to STDOUT

import argparse

def cov_from_bed(bed_file_path):

    print('contig\tcoverage')
    
    with open(bed_file_path, 'r') as bed:
        lines = bed.readlines()
        
        contig_now = lines[0].split('\t')[0]
        total_coverage = 0
        
        for line in lines:
            contig = line.split('\t')[0]
            
            if contig == contig_now:
                coverage = int(line.split('\t')[1])
                base = int(line.split('\t')[2])
                length = int(line.split('\t')[3])
                
                total_coverage += coverage*base
            else:
                cov_out = total_coverage/length
                print('{}\t{}'.format(contig_now, cov_out))
                
                contig_now = contig
                
                coverage = int(line.split('\t')[1])
                base = int(line.split('\t')[2])
                length = int(line.split('\t')[3])
                
                total_coverage = coverage*base
                
            # last fragment
            # issue: the last fragment in bed file is genome, what is this??
            if line == lines[-1]:
                cov_out = total_coverage/length
                print('{}\t{}'.format(contig, cov_out))



def main(args):

    cov_from_bed(args.input)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-b', '--input',
                        help='input the bed file',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)
