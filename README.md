# CPStools Usage

Please install Biopython packages before use.

```python
pip install biopython
```



## 1 gbcheck function

### Parameters

gbcheck.py accepts two parameters:

​	-i, --test_file

​	-r, --ref_file

The parameter of -i/--test_file is required, and -r/--ref_file is optional.

If only -i/--TEST_FILE, the gbcheck.py will check the annotation of the input file. The script will check the gene name in each gene, CDS, tRNA and rRNA. For CDS, tRNA and rRNA, it will also check the 'product' label. 

For CDS, the codon amino acids will be check, if start codon, stop codon and internal codon was wrong, it will print out the location of CDS.

### Usage

```python
python gbcheck.py -i input_file
```

The input_file is genbank format file.  And should end with '.gb'

it will check the simple wrong in the input_file.

```Python
python gbcheck.py -i input_file -r ref_file
```

The input_file and ref_file are genbank format file.  And should end with '.gb'

it will compare gene counts and difference in two GenBank files.



## 2 Sequence adjustment function

There are two scripts in this function:

​	IR.py : identify the four regions in chloroplast genome sequences.

​	seq_adj.py: this script has three modes, LSC, SSC, and RP

​		LSC: adjust the sequence start with the first bp in LSC region.

​		SSC: adjust the direction of SSC region.

​		RP: adjust the sequence to the reversed complement.

### Parameters

#### IR.py

IR.py accepts one parameter:

​	-i, --input_file

The parameter of -i/--input_file is required. The input_file can be Fasta or Genbank format file.

#### seq_adj.py

This scripts should combine the results of collinearity, which can get from the tool of nucmer.

seq_adj.py accepts four parameters, and all parameters are all required:

​	-i, --work_dir

​	-o, --save_dir

​	-f, --info_file

​	-m, --mode 

The parameter of -i/--work_dir is the directory of  sequences, which need to reverse the SSC region.

The parameter of -o/--save_dir is the directory of  sequences saved.

The parameter of -f/--info_file is the results, get from IR.py and adjusted to one line for each sequences, and separate with tab.

The parameter of -m/--mode is the adjusted mode, and should be selected form the SSC, LSC or RP. For SSC, it will adjust the direction of SSC forward, LSC to adjust the sequence start to first bp in LSC region, and RP to reverse complement the sequence and adjust to the first bp in LSC region.

### Usage

#### IR.py

```Python
python IR.py -i input_file 
```

The input_file can be Genbank of Fasta format file. it will find the four regions in the chloroplast sequences, which the IR region should surpass 1,000 bp.

#### seq_adj.py

```
python LSC_adj.py -i work_dir -o save_dir -f info.txt -m mode
```

The adjusted sequences will be save into save_dir.  -i must be directory -m must be LSC, SSC or RP.



## 3 Pi analysis function

There are two scripts in this function:

​	Pi_1.py : Extract common intergenic and gene sequences from Genbank files.

​	Pi_3.py : Calculate Pi values and sort as cp genome order.

### Parameters

#### Pi_1.py

Pi_1.py accepts one parameter:

​	-i, --work_dir: the directory of Genbank files, which need to extract IGS sequences.

Note : the Genbank files should checked by gbcheck.py

#### Pi_2.py

Pi_2.py accepts three parameter:

-i,--input: Input the directory path of the multi-alignment sequences

-r, --reference Input the file path of 'cp_sort_IGS.txt', which can generated from Pi_1.py.

-m, --mode  it has two modes, gene and IGS. IGS for sorting intergenic spaces sequences into cp order; gene for sorting gene sequences into cp order;

### Usage

#### Pi_1.py

```python
python Pi_1.py -i genbank_directory
```

The input is a directory, which contains all the genbank need for Pi analysis. This script will generate two kinds of file, 'txt' and '_IGS.fasta'. In the txt file, it records the location of intergenic, while in IGS.fasta, it is the fasta sequences of IGS. 

Then move the IGS_fasta into one directory

```shell
mkdir fasta
mv *IGS.fasta ./fasta
```

before run Pi_2.py:

Please run the command to do the multiple sequence alignment:

```shell
# IGS multi-alignment
cd unalign_common_IGS
mkdir align_IGS
for i in ./*.fasta; do mafft --auto  $i > ./align_IGS/$i ;done

# gene multi-alignment
cd common_gene/
mkdir align_gene
for i in ./*.fasta; do mafft --auto  $i > ./align_gene/$i ;done'
```

#### Pi_2.py

To calculate the Pi values and sort into the order as cp genome.

```shell
python Pi_2.py -i align_IGS -r cp_sort_IGS.txt -m gene/IGS
```

Note: if the 'rps12' was in the first line, you can move it into the correct location.

## 4 RSCU function

There are one script in this function:

​	RSCU.py: To get RSCU values from genbank files

This script will run as the follow criterion:

​	1: delete the CDS not start with 'ATG';

​	2: Only take one duplicate gene, and delete the pseudo gene

​	3: The CDS length less than --filter_length will also be delete, and filter_length was 

​		set default is 300 bp

Note: before run RSCU.py, the genbank files should be checked by gbcheck.py

### Parameters

RSCU.py accepts two parameters:

​	-i/--work_dir :  Directory of genbank files

​	-l/--filter_length: CDS filter length, default is 300 bp

### Usage

```shell
# Use default filter_length, 300 bp
python RSCU.py -i genbank_directory
# Use specified filter_length, number
python RSCU.py -i genbank_directory -l number
```

For each genbank file, it will generate a directory in the same path of genbank directory. And each genbank file will touch six file.

​	1. remove_duplicated.fasta	:	delete one duplicated gene

​	2. filter_sequences.txt	: 	the sequences was filter out.

​	3. after_filter_300.fasta	:	the saved seuquences surpass filter_length 

​	4. save_sequences_name.txt	:	the saved CDS name

​	5. merge.fasta	:	merge all the sequences into one sequence.	

​	6. RSCU_results.txt:	calculate RSCU values (table=11)

## 5 SSRs function

There are one script in this function:

​	SSR_analysis.py: To find SSRs  from fasta/genbank files

Note: to ensure the accuracy of the results, we should perform collinearity analysis on all sequences, which can be completed under the sequence adjustment function

### Parameters

SSR_analysis.py accepts two parameters:

​	-i/--input_file :  the input file can be Genbank or Fasta format

​	-k/--kmer_lengths: the length of the six types of SSR, default is 10,6,5,4,4,4

### Usage

```shell
# Use default kmer_lengths, 10,6,5,4,4,4
python SSR_analysis.py -i input_file
# Use specified kmer_lengths, number 11,6,5,4,4,3
python SSR_analysis.py -i input_file -k 11,6,5,4,4,3
```



## 6 Phylogenetic analysis function

There are two mode in this function:

-m/--mode cds/pro cds for common cds sequences; pro for common
                    protein sequences

### Parameters

common_cds.py and common_pro.py both accept two parameters:

​	-i/--ref_file: All the Genbank files was set at a directory, and choose one file as ref_file

​	-o/--output_dir: output directory of fasta files

​	-m/--mode cds/pro  cds for common cds sequences; pro for common
​                protein sequences

### Usage

```shell
# The merged cds file is saved in results/merge_cds.fasta
Python phy_built.py -i ref_file.gb -o results -m cds
# The merged protein file is saved in results2/merge_pro.fasta
Python phy_built.py -i ref_file.gb -o results2 -m pro
```



## 7 Format conversion function

There are three modes in this function:

​	gb2fa: To convert Genbank format file into Fasta format.

​	gb2mVISTA: To prepare annotation file for mVISTA software.

​	gb2tbl: To convert Genbank format file into tabular format.

Note: to ensure the accuracy of the results, the Genbank files should checked by gbcheck.py

### Parameters

accepts three parameters:

​	-i/--gb_file:  Input path of genbank file

​	-o/--fa_file: output path of fasta file

​	-m/--mode: fasta,mVISTA,tbl

### Usage

#### gb2fa

```shell
python converse.py -i input.gb -o output.fasta -m fasta
```

#### gb2tbl.py

```shell
# the genbank is the directory of Genbank format files
python converse.py -i genbank -o output.fasta -m tbl
```

#### gb2mVISTA.py

```shell
# the genbank is the directory of Genbank format files
python converse.py -i genbank -o output.fasta -m mVISTA
```















































