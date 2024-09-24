from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import pickle
from Bio import SeqIO
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='This is a script for processing files.')
    parser.add_argument('-i', '--input', type=str, default='examples/examples.fasta')
    parser.add_argument('-o', '--output', type=str, default='result/result.csv')
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    return input_file, output_file

def samples(samples_path):
	samples_seq={}
	for record in SeqIO.parse(samples_path, "fasta"):
		seq = str(record.seq)
		gene_name = record.description
		if gene_name not in samples_seq:
			samples_seq[gene_name] = [seq]
		else:
			samples_seq[gene_name].append(seq)
	print(f'Samples upload finished! Gene number:{len(samples_seq)}')
	return samples_seq

def feature_seq(seq):
    '''Generate k-mer features of the sequence'''
    kmax=3
    λmax=1
    list_feature=[]
    kk = range(0,kmax)
    λλ = range(0,λmax+1)
    for λ in λλ:
        for k in kk:
            sequences_1 = []
            sequences_2 = []
            sequences_3 = []
            sequences = []
            n=0
            short_seqs = ['']
            base = ['A','T','C','G']
            while True:
                if n >= (4**(k+1)):
                    break
                else:
                    short_seq_temp = []
                    for b in base:
                        for short_seq in short_seqs:
                            short_seq = b + short_seq
                            short_seq_temp.append(short_seq)
                            n += 1
                    short_seqs = short_seq_temp
            for j in range(0,len(seq),3):
                i_1 = []
                i_2 = []
                i_3 = []
                if k == 0:
                    if λ == 0:
                        for i in range(0,k+1): 
                            index_1 = 0+i+j
                            index_2 = 1+i+j
                            index_3 = 2+i+j
                            i_1.append(index_1)
                            i_2.append(index_2)
                            i_3.append(index_3)
                else: 
                    for i in range(0,k):
                        index_1 = 0+i+j
                        index_2 = 1+i+j
                        index_3 = 2+i+j
                        i_1.append(index_1)
                        i_2.append(index_2)
                        i_3.append(index_3)
                    index_1 = 0+k+λ+j
                    index_2 = 1+k+λ+j
                    index_3 = 2+k+λ+j
                    i_1.append(index_1)
                    i_2.append(index_2)
                    i_3.append(index_3)
                seq_1=''
                seq_2=''
                seq_3=''
                for i in i_1:
                    if i in range(0,len(seq)):
                        seq_1 += seq[i]
                        if len(seq_1) == k+1:
                            sequences_1.append(seq_1)
                for i in i_2:
                    if i in range(0,len(seq)):
                        seq_2 += seq[i]
                        if len(seq_2) == k+1:
                            sequences_2.append(seq_2)
                for i in i_3:
                    if i in range(0,len(seq)):
                        seq_3 += seq[i]
                        if len(seq_3) == k+1:
                            sequences_3.append(seq_3)
            sequences.append(sequences_1)
            sequences.append(sequences_2)
            sequences.append(sequences_3)
            for sequence in sequences:
                if sequence: 
                    for short_seq in short_seqs:
                        ratio = sequence.count(short_seq)/len(sequence)
                        ratio =float ('%.5f' % ratio)
                        list_feature.append(ratio)
    return list_feature

def feature(samples_seq):
    df_E = pickle.load(open('model/Expr_Feature.pkl','rb'))
    #Generate sequence features
    gene_lists = []
    features = {}
    for k,v in samples_seq.items():
        feature = []
        for seq in v:
            feature.append(feature_seq(seq))
        df_temp = pd.DataFrame(feature)
        column_means = df_temp.mean()
        features[k]=column_means
    df_seq = pd.DataFrame(features).T
    feature_df = pd.merge(df_E, df_seq, left_index=True, right_index=True)
    print(f'\nFeature matrix generation completed! Matrix shape:{feature_df.shape}')
    missing_genes = list(set(df_seq.index.tolist())-set(feature_df.index.tolist()))
    if missing_genes:
        print(f'missing_genes:{",".join(missing_genes)}. These genes are not present in the gene expression table, therefore we cannot generate their corresponding features.')
    return feature_df

def main(input_file,output_file):
    #Load samples
    samples_path = input_file
    samples_seq = samples(samples_path)
    #Feature generation
    feature_df = feature(samples_seq)
    selector = pickle.load(open('model/selector.pkl','rb'))
    feature_df= feature_df.iloc[:,selector]
    print(feature_df.shape)
    feature_df.columns = feature_df.columns.astype(str)
    #prediction
    model = pickle.load(open('model/rf.pkl','rb'))
    prediction = model.predict(feature_df)
    probability= model.predict_proba(feature_df)[:, 1]
    df_result = pd.DataFrame({'prediction':prediction,'predict_proba':probability},index=feature_df.index.tolist())
    df_result.to_csv(output_file)
    print(f'\nFinished\n')

if __name__ == '__main__':
    input_file, output_file = parse_arguments()
    main(input_file,output_file)
