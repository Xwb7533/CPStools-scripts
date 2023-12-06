from Bio import SeqIO
import os
import sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Get mVISTA information from genbank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def gb2mVISTA(input_file):


    all_info = []
    for rec in SeqIO.parse(input_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'gene':
                for part in feature.location.parts:
                    if int(part.strand) == 1:
                        if ['>', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]] not in all_info:
                            all_info.append(['>', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]])
                    else:
                        if ['<', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]] not in all_info:
                            all_info.append(['<', int(part.start) + 1, int(part.end), feature.qualifiers['gene'][0]])
            elif feature.type == 'CDS':
                for part in feature.location.parts:
                    if int(part.strand) == 1:
                        if [int(part.start) + 1, int(part.end), 'exon'] not in all_info:
                            all_info.append([int(part.start) + 1, int(part.end), 'exon'])
                    else:
                        if [int(part.start) + 1, int(part.end), 'exon'] not in all_info:
                            all_info.append([int(part.start) + 1, int(part.end), 'exon'])
            elif feature.type == 'rRNA' or feature.type == 'tRNA':
                for part in feature.location.parts:
                    if int(part.strand) == 1:
                        all_info.append([int(part.start) + 1, int(part.end), 'utr'])
                    else:
                        all_info.append([int(part.start) + 1, int(part.end), 'utr'])
            else:
                pass
    # write results into save_file
    dir_name = os.path.dirname(input_file)
    file_name = os.path.basename(input_file).split('.')[0]
    save_file = os.path.join(dir_name, file_name + '_mVISTA.txt')
    with open(save_file, 'w') as ff:
        for list_ in all_info:
            if len(list_) == 4:
                ff.write(f"{list_[0]} {list_[1]} {list_[2]} {list_[3]}\n")
            else:
                ff.write(f"{list_[0]} {list_[1]} {list_[2]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of genbank files")
    args = parser.parse_args()
    if args.work_dir:
        for i in os.listdir(args.work_dir):
            if i.endswith('gb') or i.endswith('gbk'):
                file_path = os.path.join(args.work_dir, i)
                print(file_path)
                gb2mVISTA(file_path)
            else:
                print(f"{i} is not genbank format files, Please endswith 'gb' or 'gbk'")
        print("Format converse ok!")
    else:
        parser.print_help()


