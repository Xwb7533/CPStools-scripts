import os
import argparse
from Bio import SeqIO


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Calculate Pi values
        The results is same as DNAsp6!
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def calculate_Pi_values(work_dir):
    all_pi_results = []
    for align_fasta in os.listdir(work_dir):
        if align_fasta.endswith('.fasta'):
            a, b, c, d = [], [], [], []
            pi = 0
            fasta_file = os.path.join(work_dir, align_fasta)
            for rec in SeqIO.parse(fasta_file, format='fasta'):
                a.append(rec.id)
                for i in range(len(rec.seq)):
                    if rec.seq[i] == '-':
                        b.append(i)
            all_number = len(a) * (len(a) - 1) / 2
            # delete all have '-' in seq location
            all_del = sorted(set(b))[::-1]
            for rec in SeqIO.parse(fasta_file, 'fasta'):
                for x in all_del:
                    rec.seq = list(rec.seq)
                    del rec.seq[x]
                d.append(rec.seq)
            # statistics same and diff
            for y in range(len(d[0])):
                c = []
                for x in d:
                    c.append(x[y])
                diff = 0
                for sig in range(len(a)):
                    for sig2 in range(sig + 1, len(a)):
                        if c[sig] != c[sig2]:
                            diff += 1
                pi += diff / all_number
            if len(d[0]) == 0:
                final_pi = 0
            else:
                final_pi = format(pi / len(d[0]), '.5f')
            all_pi_results.append(f"{align_fasta[:-6]}\t{final_pi}")
            print(f"{align_fasta[:-6]}\t{final_pi}")
    with open(os.path.join(work_dir, 'Pi_results.txt'), 'w') as ff:
        for each_pi in all_pi_results:
            ff.write(f"{each_pi}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input_dir', type=str, help='Align fasta format work directory')
    args = parser.parse_args()
    if args.input_dir:
        work_dir = os.path.abspath(args.input_dir)
        if not os.path.isdir(work_dir):
            print("Please input align fasta format work directory, not fasta file")
        else:
            calculate_Pi_values(work_dir)
            Pi_results = os.path.abspath(os.path.join(work_dir, 'Pi_results.txt'))
            print(f"The calculate results have been written into: \n {Pi_results} \n")
    else:
        parser.print_help()