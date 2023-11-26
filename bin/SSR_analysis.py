import re
import os
import argparse
from Bio import SeqIO


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Find SSRs in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com
        """

        return f"{description}\n\n{help_text}"


def parse_file(input_file):
    for rec in SeqIO.parse(input_file, 'fasta'):
        return str(rec.seq)


def find_SSRs(input_file, *length):
    type_length = []
    if len(length) == 0:
        type_length = [10, 6, 5, 4, 4, 4]
    elif len(length) == 6:
        type_length = length
    else:
        print("Please input type length equal 6 or 0")
        return []
    k1, k2, k3, k4, k5, k6 = type_length
    print(k1, k2, k3, k4, k5, k6)
    print(f"The parameter is set as:\n"
          f"Mononucleotide:{k1}\n"
          f"Dinucleotide:{k2}\n"
          f"Trinucleotide:{k3}\n"
          f"Tetranucleotide:{k4}\n"
          f"Pentanucleotide:{k5}\n"
          f"Hexanucleotide:{k6}\n")

    my_seq = parse_file(input_file)
    match1 = re.finditer(r'([ATCG])\1{{{0},}}'.format(k1-1), my_seq)
    match2 = re.finditer(r'(((?!AA|TT|CC|GG)[ATCG]){{2}})\1{{{0},}}'.format(k2-1), my_seq)
    match3 = re.finditer(r'(((?!AAA|TTT|CCC|GGG)[ATCG]){{3}})\1{{{0},}}'.format(k3-1), my_seq)
    match4 = re.finditer(r'(?!(\w)(\w)\1\2)(((?!AAAA|TTTT|CCCC|GGGG)[ATCG]){{4}})\3{{{0},}}'.format(k4-1), my_seq)
    #match4 = re.finditer(r'(((?!AAAA|TTTT|CCCC|GGGG)[ATCG]){{4}})\1{{{0},}}'.format(k4-1), my_seq)
    match5 = re.finditer(r'(((?!AAAAA|TTTTT|CCCCC|GGGGG)[ATCG]){{5}})\1{{{0},}}'.format(k5-1), my_seq)
    match6 = re.finditer(r'(?!(\w)(\w)(\w)\1\2\3)(((?!AAAAAA|TTTTTT|CCCCCC|GGGGGG)[ATCG]){{6}})\4{{{0},}}'.format(k6-1), my_seq)
    #match6 = re.finditer(r'(((?!AAAAAA|TTTTTT|CCCCCC|GGGGGG)[ATCG]){{6}})\1{{{0},}}'.format(k6-1), my_seq)

    print(k1-1, k2-1, k3-1, k4-1, k5-1, k6-1)
    # 收集所有匹配结果
    all_matches = []
    for match in match1, match2, match3, match4, match5, match6:
        all_matches.extend(match)
    
    print(f"results:\ttotal {len(all_matches)} were detected! and the results was written into SSRs_results.txt")
    # 按照span进行排序
    all_matches = sorted(all_matches, key=lambda m: m.span())
    cc, dd = [], []
    for i in all_matches:
        cc.append(i[0])
    file_name = os.path.basename(input_file).split('.')[0]
    save_name = file_name + "_SSRs_results.txt"
    save_name2 = file_name + "_SSRs_simplification.txt"
    file_save = os.path.join(os.path.dirname(input_file), save_name)
    file_save2 = os.path.join(os.path.dirname(input_file), save_name2)
    with open(file_save, 'w') as ff:
        ff.write("type\tlength\tstart\tend\n")
        for match_obj in all_matches:
            repeat_unit = match_obj.group()  # 获取匹配到的具体内容
            start_pos = match_obj.start()
            end_pos = match_obj.end()
            ssr_type = len(repeat_unit)  # SSR类型为匹配到的具体内容的长度
            ff.write(f"{repeat_unit}\t{len(repeat_unit)}\t{start_pos+1}\t{end_pos}\n")

def main():
    # print("start")
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--fasta_file', help='fasta format file', required=True)
    parser.add_argument('-k', '--kmer_length', help='SSRs length, default is 10,6,5,4,4,4')
    args = parser.parse_args()
    input_file = args.fasta_file
    kmer_lengths = list(map(int, args.kmer_length.split(','))) if args.kmer_length else []
    result = find_SSRs(input_file, *kmer_lengths)
        


if __name__ == "__main__":
    main()