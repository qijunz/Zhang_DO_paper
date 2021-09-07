# Program to print a table of sequences in a multifasta file and their gc content
# USAGE: fasta_gc_table.py -s <input fasta> [> out]
# prints to STDOUT

import argparse

def length_calculation(fasta):

    print('Sequence_name\tlength')
    
    with open(fasta, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith('>'):
                if not line == lines[0]:
                    print('{}\t{}'.format(seq_name, length))
                seq_name = line.strip().split('>')[-1]
                length = 0
            else:
                length += int(len(line.strip()))
        print('{}\t{}'.format(seq_name, length))


def main(args):

    length_calculation(args.input)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-s', '--input',
                        help='input the fasta file',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)

