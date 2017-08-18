from pandas import DataFrame as df
import csv
def construct_specific_snv(test_csv,output_snv):
    info_df = df.from_csv(test_csv,index_col=False)
    count = 0
    gene_list = info_df.loc[:,'Gene.refGene'].values.tolist()
    gene_list = [_i.split(';') for _i in gene_list]
    with open(output_snv,'w') as f1:
        spamwriter = csv.writer(f1)
        for ks in gene_list:
            for k in ks:
                spamwriter.writerow([count,k])
            count +=1