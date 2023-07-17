# Prediction of Gene Essentiality in Normal Human Individuals

#### Description
The model utilizes both gene expression and DNA sequence features to predict the essentiality of human genes using a random-forest approach.

#### Requirements
1.python == 3.8
2.sklearn == 0.24.1
3.pandas == 1.2.4
4.Bio == 1.79

#### How to run our code
1.Please provide a fasta format file of your DNA sequence, with the gene symbol after '>' in the file.
eg.
>PF4
ATGAGCTCCGCAGCCGGGTTCTGCGCCTCACGCCCCGGGCTGCTGTTCCTGGGGTTGCTGCTCCTGCCACTTGTGGTCGC
CTTCGCCAGCGCTGAAGCTGAAGAAGATGGGGACCTGCAGTGCCTGTGTGTGAAGACCACCTCCCAGGTCCGTCCCAGGC
ACATCACCAGCCTGGAGGTGATCAAGGCCGGACCCCACTGCCCCACTGCCCAACTGATAGCCACGCTGAAGAATGGAAGG
AAAATTTGCTTGGACCTGCAAGCCCCGCTGTACAAGAAAATAATTAAGAAACTTTTGGAGAGTTAG
2.Download the gene expression matrix from (https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz) and put it into '/data' folder.
3.Change the file paths of input file and output file.
4.Execute main.py to run the program, and the results will be outputted in the '/result' folder.
