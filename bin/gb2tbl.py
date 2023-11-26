from Bio import SeqIO
import os
import sys
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Get tbl information from genbank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


output_order = ["trans_splicing", "exception", "pesudo", "codon_start", "product", "gene", "transl_table"]


def gb2tbl(input_file):
    dir_name = os.path.dirname(input_file)
    file_name, file_extension = os.path.splitext(os.path.basename(file_path))
    save_file = os.path.join(dir_name, file_name + '.tbl')
    with open(save_file, 'w') as ff:
        for rec in SeqIO.parse(input_file, 'gb'):
            ff.write(f">Feature {rec.id}\n")
            for feature in rec.features[1:]:
                if len(feature.location.parts) == 1:
                    for part in feature.location.parts:
                        if part.strand == 1:
                            ff.write(f"{part.start + 1}\t{part.end}\t{feature.type}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                        else:
                            ff.write(f"{part.end}\t{part.start + 1}\t{feature.type}\n")
                            for key in output_order:
                                if key not in feature.qualifiers:
                                    continue
                                value = feature.qualifiers[key][0]
                                ff.write(f"\t\t\t{key}\t{value}\n")
                elif len(feature.location.parts) == 2:
                    if feature.location.strand == -1:
                        ff.write(f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n"
                              f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
                    elif feature.location.strand == 1:
                        ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                              f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
                    else:
                        ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                              if feature.location.parts[0].strand == 1 else
                              f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n")
                        ff.write(f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                              if feature.location.parts[1].strand == 1 else
                              f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
                else:
                    if feature.location.strand == -1:
                        ff.write(f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n"
                              f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n"
                              f"{feature.location.parts[2].end}\t{feature.location.parts[2].start + 1}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
                    elif feature.location.strand == 1:
                        ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                              f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                              f"{feature.location.parts[2].start + 1}\t{feature.location.parts[2].end}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
                    else:
                        ff.write(f"{feature.location.parts[0].start + 1}\t{feature.location.parts[0].end}\t{feature.type}\n"
                              if feature.location.parts[0].strand == 1 else
                              f"{feature.location.parts[0].end}\t{feature.location.parts[0].start + 1}\t{feature.type}\n")
                        ff.write(f"{feature.location.parts[1].start + 1}\t{feature.location.parts[1].end}\n"
                              if feature.location.parts[1].strand == 1 else
                              f"{feature.location.parts[1].end}\t{feature.location.parts[1].start + 1}\n")
                        ff.write(f"{feature.location.parts[2].start + 1}\t{feature.location.parts[2].end}\n"
                              if feature.location.parts[2].strand == 1 else
                              f"{feature.location.parts[2].end}\t{feature.location.parts[2].start + 1}\n")
                        for key in output_order:
                            if key not in feature.qualifiers:
                                continue
                            value = feature.qualifiers[key][0]
                            ff.write(f"\t\t\t{key}\t{value}\n")
        ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of genbank files")
    args = parser.parse_args()
    if args.work_dir:
        for i in os.listdir(args.work_dir):
            if i.endswith('gb') or i.endswith('gbk'):
                file_path = os.path.abspath(os.path.join(args.work_dir, i))
                file_name2, file_extension = os.path.splitext(os.path.basename(file_path))
                save_file2 = os.path.abspath(os.path.join(args.work_dir, file_name2 + '.tbl'))
                print(f"The input file is:\n\t\t{file_path}\nand the tbl results have been saved in:\n"
                    f"\t\t{save_file2}\n{'-' * 40}\n")
                gb2tbl(file_path)
            else:
                print(f"{i} is not genbank format files, Please endswith 'gb' or 'gbk'\n{'-' * 40}\n")
        print("Format converse ok!")
    else:
        parser.print_help()
