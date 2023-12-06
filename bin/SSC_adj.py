import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """Adjust the LSC start in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def adjust_SSC_forward(work_dir, save_dir, info_text):
    info_cont = open(info_text, 'r')
    info_line = info_cont.readline()

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    while info_line:
        file_name = info_line.split('\t')[0]
        SSC_loc_start = int(info_line.split('\t')[3].split('-')[0].split(':')[1]) - 1
        SSC_loc_end = int(info_line.split('\t')[3].split('-')[1].strip())
        print(file_name)
        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, file_name)
        if os.path.exists(work_file):
            with open(save_file, 'w') as ff:
                for rec in SeqIO.parse(work_file, 'fasta'):
                    SSC_seq = Seq(str(rec.seq[SSC_loc_start:SSC_loc_end])).reverse_complement()
                    rev_seq = str(rec.seq[:SSC_loc_start]) + str(SSC_seq) + str(rec.seq[SSC_loc_end:])
                    ff.write(f'>{rec.id}\n{rev_seq}\n')
                info_line = info_cont.readline()
        else:
            abs_path = os.path.abspath(work_file)
            print(f"No such file: {abs_path}")
            info_line = info_cont.readline()
    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of fasta files")
    parser.add_argument("-o", "--save_dir", help="Output directory for save files")
    parser.add_argument("-f", "--info_file", help="file of the information of the fasta file four regions")
    args = parser.parse_args()
    if args.work_dir and args.save_dir and args.info_text:
        adjust_SSC_forward(args.work_dir, args.save_dir, args.info_file)
    else:
        parser.print_help()