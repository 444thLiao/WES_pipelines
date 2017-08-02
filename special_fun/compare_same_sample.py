from pandas import DataFrame as df
from matplotlib_venn import venn2_unweighted

def gene_filter(file_df,_in = []):
    #from filters.py
    gene_var_bucket = []
    indexs = file_df.index
    for _idx in list(indexs):
        if file_df.loc[_idx,'Gene.refGene']:
            genes = file_df.loc[_idx,'Gene.refGene'].split(';')
            for _g in genes:
                if _g in _in and _idx not in gene_var_bucket:
                    gene_var_bucket.append(_idx)
    return gene_var_bucket

with open('/home/liaoth/project/YZ_XR/genes_355.list') as f1:
    gene_list = f1.read().split('\n')

def init_compare(input1,input2,labels=None):
    input1_df = df.from_csv(input1, index_col=False)
    input2_df = df.from_csv(input2, index_col=False)
    input1_df_index = [';'.join(_i.split('\t')[:5]) for _i in list(input1_df.Otherinfo)]
    input2_df_index = [';'.join(_i.split('\t')[:5]) for _i in list(input2_df.Otherinfo)]
    if labels:
        ax = venn2_unweighted([set(input1_df_index),set(input2_df_index)],set_labels=labels)
        return ax
    else:
        input1_df.index = input1_df_index
        input2_df.index = input2_df_index
        result = {}

        result['small_one_unique'] = input1_df.loc[set(input1_df_index).difference(set(input2_df_index)),:]
        result['big_one_should_in_small'] = input2_df.loc[gene_filter(input2_df,_in=gene_list),:]
        result['shared'] = input1_df.loc[set(input1_df_index).intersection(set(input2_df_index)),:]
        return result


stats={}
stats_index = {}
stats['Total tumor coverage = 0'] = len(a[a.loc[:,['T_mut_cov','N_ref_cov']].sum(1)==0])

stats_index['N_mut_per > 3%'] = list(a[a.loc[:,'N_mut_per'] > 0.03].index)
stats_index['Total tumor coverage < 10'] = list(a[a.loc[:,['T_mut_cov','T_ref_cov']].sum(1)<10].index)
stats_index['Total tumor coverage = 0'] = list(a[a.loc[:,'T_mut_per']==0].index)
stats_index['N_mut_per = T_mut_per = 0'] = list(a[(a.loc[:,'N_mut_per'] == a.loc[:,'T_mut_per']) & (a.loc[:,'T_mut_per'] == 0)].index)
stats_index['N_mut_per = T_mut_per = 1'] = list(a[(a.loc[:,'N_mut_per'] == a.loc[:,'T_mut_per']) & (a.loc[:,'T_mut_per'] == 1)].index)
stats_index['N_mut_per > T_mut_per'] = list(a[a.loc[:,'N_mut_per'] > a.loc[:,'T_mut_per']].index)
