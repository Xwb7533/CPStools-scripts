import os
from Bio import SeqIO
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Intergenic Sequences Extraction.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def IGS_extract(input_file):
    for rec in SeqIO.parse(input_file, format='genbank'):
        genome_length = [[int(part.end)] for part in rec.features[0].location.parts][0][0]
        my_seq = rec.seq
        all_feature = []
        all_info = []
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                all_feature.append(feature)
        for i in range(len(all_feature)):
            #print(input_file)
            # print(all_feature[i])
            gene_name = all_feature[i].qualifiers['gene'][0]
            gene_location = all_feature[i].location.parts
            gene1_exon_info = [[int(part.start), int(part.end), part.strand] for part in gene_location]
            exon_number = len(gene_location)
            if exon_number == 1:
                all_info.append(f"{gene_name}\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
            if exon_number == 2:
                if gene1_exon_info[0][1] == genome_length:
                    all_info.append(f"{gene_name}\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[0][2]}")
                else:
                    all_info.append(f"{gene_name}_1\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
                    all_info.append(f"{gene_name}_2\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[1][2]}")
            if exon_number == 3:
                all_info.append(f"{gene_name}_1\t{gene1_exon_info[0][0]}\t{gene1_exon_info[0][1]}\t{gene1_exon_info[0][2]}")
                all_info.append(f"{gene_name}_2\t{gene1_exon_info[1][0]}\t{gene1_exon_info[1][1]}\t{gene1_exon_info[1][2]}")
                all_info.append(f"{gene_name}_3\t{gene1_exon_info[2][0]}\t{gene1_exon_info[2][1]}\t{gene1_exon_info[2][2]}")
        all_info = list(set(all_info))
        all_info.sort(key=lambda x: int(x.split('\t')[1]))
        save_file = os.path.join(os.path.abspath(input_file).split('.')[0] + '_intergenic_location.txt')
        save_file_w = open(save_file, 'w')
        for i in range(len(all_info)-1):
            info_list = all_info[i].split('\t')
            next_list = all_info[i+1].split('\t')
            save_file_w.write(f"{info_list[0]}-{next_list[0]}\t{info_list[-2]}\t{next_list[1]}\n")
        end_gene_info = all_info[-1].split('\t')
        start_gene_info = all_info[0].split('\t')
        if int(end_gene_info[-2]) < int(start_gene_info[1]):
            # print(f"{end_gene_info[0]}-{start_gene_info}[0]\t{end_gene_info[-2]}\t{start_gene_info[1]}")
            save_file_w.write(f"{end_gene_info[0]}-{start_gene_info}[0]\t{end_gene_info[-2]}\t{start_gene_info[1]}\n")
        else:
            if int(end_gene_info[2]) < genome_length:
                # print(f"{end_gene_info[0]}-{start_gene_info[0]}\t{end_gene_info[-2]}\t{genome_length}\t0\t\
                # {start_gene_info[1]}")
                save_file_w.write(f"{end_gene_info[0]}-{start_gene_info[0]}\t{end_gene_info[-2]}\t{genome_length}\t0\t\
                {start_gene_info[1]}\n")
            else:
                pass
        save_file_w.close()
        all_fasta_file = os.path.join(os.path.abspath(input_file).split('.')[0] + '_IGS.fasta')
        all_fasta = open(all_fasta_file, 'w')
        save_results = open(save_file, 'r')
        result_line = save_results.readline().strip()
        while result_line:
            result_line_list = result_line.split('\t')
            if len(result_line_list) == 3:
                if int(result_line_list[2]) > int(result_line_list[1]):
                    all_fasta.write(f">{result_line_list[0]}\n{my_seq[int(result_line_list[1]):int(result_line_list[2])]}\n")
                    result_line = save_results.readline().strip()
                else:
                    print(f"{result_line_list[0]} has overlap!")
                    result_line = save_results.readline().strip()

            else:
                all_fasta.write(f">{result_line_list[0]}\n{my_seq[int(result_line_list[1]):int(result_line_list[2])]}\
                {my_seq[int(result_line_list[3]):int(result_line_list[4])]}\n")
                result_line = save_results.readline().strip()
        all_fasta.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", required=True, help="Input directory of genbank files")
    args = parser.parse_args()
    if not args.work_dir:
        parser.error("Please provide the input directory of genbank files using -i or --work_dir")
    else:
        for i in os.listdir(args.work_dir):
            if i.endswith('gb') or i.endswith('gbk'):
                file_path = os.path.join(args.work_dir, i)
                print(file_path)
                IGS_extract(file_path)
            else:
                print(f"please input genbank format files and endwith 'gb' or 'gbk'! ")

print("All Done !!!")
