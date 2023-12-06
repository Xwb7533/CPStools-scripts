from Bio import SeqIO
import os
import sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        translate genbank format files to fasta format.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def parse_genbank_file(input_file):
    try:
        with open(input_file, 'r') as f:
            records = list(SeqIO.parse(f, 'genbank'))
            if records:
                return records
            else:
                print(f"No valid GenBank records found in {input_file}")
                return None
    except (ValueError, IndexError, KeyError, IOError, FileNotFoundError) as e:
        print(f"Failed to parse GenBank file '{input_file}': {e}")
        return None


def gb2fa(input_file, save_file):
    if parse_genbank_file(input_file):
        with open(save_file, 'w') as ff:
            for rec in SeqIO.parse(input_file, format='genbank'):
                ff.write(f'>{rec.id}\n{rec.seq}\n')
            ff.close()
    else:
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--gb_file", help="Input path of genbank file")
    parser.add_argument("-o", "--fa_file", help="output path of fasta file")
    args = parser.parse_args()
    if args.gb_file and args.fa_file:
        if args.gb_file.endswith('gb') or args.gb_file.endswith('gbk'):
            file_path = os.path.abspath(args.gb_file)
            save_file = os.path.abspath(args.fa_file)
            print(file_path)
            gb2fa(file_path, save_file)
            print("Format converse ok!")
        else:
            print(f"{args.gb_file} is not genbank format files, Please endswith 'gb' or 'gbk'")
    else:
        parser.print_help()
    


