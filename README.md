# Predict the intolerance of loss-of-function in human genes

#### Description
The model utilizes both gene expression and DNA sequence features to predict the intolerance of loss-of-function using a random forest-based (RF-based) approach.

#### Requirements
1.python == 3.8<br>
2.sklearn == 0.24.1<br>
3.pandas == 1.2.4<br>
4.Bio == 1.79<br>

#### How to run our code
</p>1. Please provide a fasta format file of your DNA sequence, with the gene symbol after &apos;>&apos; in the file.<br>
eg.<br>
>PF4<br>
ATGAGCTCCGCAGCCGGGTTCTGCGCCTCACGCCCCGGGCTGCTGTTCCTGGGGTTGCTGCTCCTGCCACTTGTGGTCGC
CTTCGCCAGCGCTGAAGCTGAAGAAGATGGGGACCTGCAGTGCCTGTGTGTGAAGACCACCTCCCAGGTCCGTCCCAGGC
ACATCACCAGCCTGGAGGTGATCAAGGCCGGACCCCACTGCCCCACTGCCCAACTGATAGCCACGCTGAAGAATGGAAGG
AAAATTTGCTTGGACCTGCAAGCCCCGCTGTACAAGAAAATAATTAAGAAACTTTTGGAGAGTTAG</p>
</p>2. Download the gene expression matrix from ( https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz).</p>
</p>3. Reference code:<br>
python main.py -i intput.fasta -o result.csv -g GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct<br>
-i: input file
-o: output file
-g: gene expression matirx file
</p>

