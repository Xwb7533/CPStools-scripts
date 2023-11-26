import os,sys
from Bio import SeqIO
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Extract common gene sequences from genbank files
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def common_gene_extract(input_file):
    work_dir = os.path.dirname(input_file)
    all_gene = []
    for rec in SeqIO.parse(input_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                if feature.qualifiers['gene'][0] not in all_gene:
                    all_gene.append(feature.qualifiers['gene'][0])
    print(len(all_gene))
    for files in os.listdir(work_dir):
        if files.endswith('gb') or files.endswith('gbk'):
            single_gene = []
            gb_file = os.path.join(work_dir, files)
            print(gb_file)
            for rec in SeqIO.parse(gb_file, format='genbank'):
                for feature in rec.features:
                    if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                        if feature.qualifiers['gene'][0] not in single_gene:
                            single_gene.append(feature.qualifiers['gene'][0])
                print(len(single_gene))
            # delete unigue gene
            for gene_index in range(len(all_gene)-1, -1, -1):
                if all_gene[gene_index].lower() not in [y.lower() for y in single_gene]:
                    all_gene.remove(all_gene[gene_index])
    print("Total gene number is : ", end='\t')
    print(len(all_gene))
    gene_name_file = os.path.join(os.path.dirname(work_dir), 'gene_cp_sort.txt')
    with open(gene_name_file, 'w') as ff:
        for i in all_gene:
            ff.write(f'{i}\n')
    save_dir = os.path.join(os.path.dirname(work_dir), 'common_gene')
    if os.path.exists(save_dir):
        print(
            f"\t\t########  Run failed !!!  ########\n"
            f"The save directory has been existed, please delete the directory!\n"
            f"\t\t\t{os.path.abspath(save_dir)}")
        sys.exit()
    os.mkdir(save_dir)
    for gene_name in all_gene:
        file_name = str(gene_name) + '.fasta'
        file_path = os.path.join(save_dir, file_name)
        with open(file_path, 'w') as fasta_file:
            for gb_file in os.listdir(work_dir):
                if gb_file.endswith('gb') or gb_file.endswith('gbk'):
                    gb_file_path = os.path.join(work_dir, gb_file)
                    fasta_file.write(f">{gb_file.split('.')[0]}\n")
                    for rec in SeqIO.parse(gb_file_path, format='genbank'):
                        my_seqs = []
                        for feature in rec.features:
                            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                                if feature.qualifiers['gene'][0].lower() == gene_name.lower():
                                    my_seqs.append(feature.extract(rec.seq))
                        if len(my_seqs) == 1:
                            fasta_file.write(f"{my_seqs[0]}\n")
                        if len(my_seqs) == 2:
                            my_seqs.remove(my_seqs[0]) if len(my_seqs[0]) <= len(my_seqs[1]) else my_seqs.remove(my_seqs[1])
                            fasta_file.write(f"{my_seqs[0]}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input the file path of reference genbank file')
    args = parser.parse_args()
    if not args.input_file:
        parser.error("Please provide the input file path of reference genbank file using -i or --input_file")
        sys.exit()
    common_gene_extract(os.path.abspath(args.input_file))
    save_results_dir = os.path.abspath(os.path.join(os.path.dirname(args.input_file), "common_gene"))

    print(
        f"\t\t########  Run successfully !!!  ########\n"
        f"The common genes sequences have been extracted, and the results were save into the directory:\n"
        f"\t\t\t{save_results_dir}")

    print(
    "##### Next step is to do multiple alignment, the command is : ######\n"
    "\t\t1: 'cd common_gene/ \n"
    "\t\t2: 'mkdir align_gene'\n"
    "\t\t3: 'for i in ./*.fasta; do mafft --auto  $i > ./align_gene/$i ;done'\n")
