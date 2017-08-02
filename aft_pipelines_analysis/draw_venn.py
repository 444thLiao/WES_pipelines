from matplotlib_venn import venn2_unweighted
from pandas import DataFrame as df
from matplotlib import pyplot
def venn2_draw(csv1,csv2,output,labels=[]):
    DF1 = df.from_csv(csv1,index_col=False)
    DF2 = df.from_csv(csv2,index_col=False)

    index1 = [';'.join([str(_i) for _i in list(DF1.iloc[_idx, :5])]) for _idx in list(DF1.index)]
    index2 = [';'.join([str(_i) for _i in list(DF2.iloc[_idx, :5])]) for _idx in list(DF2.index)]
    if labels:
        ax = venn2_unweighted([set(index1),set(index2)],set_labels=labels)
        ax.get_label_by_id('A').set_fontsize(25)
        ax.get_label_by_id('B').set_fontsize(25)
        for _i in ax.subset_labels:
            _i.set_fontsize(25)
        pyplot.savefig(output)


def extract_shared_snvs(csv1,csv2,output,samples = []):
    DF1 = df.from_csv(csv1,index_col=False)
    DF2 = df.from_csv(csv2,index_col=False)

    index1 = [';'.join([str(_i) for _i in list(DF1.iloc[_idx, :5])]) for _idx in list(DF1.index)]
    index2 = [';'.join([str(_i) for _i in list(DF2.iloc[_idx, :5])]) for _idx in list(DF2.index)]

    DF1.index = index1
    DF2.index = index2

    shared_index = set(index1).intersection(set(index2))

    shared_DF1 = DF1.ix[sorted(shared_index,reverse=True)]
    shared_DF2 = DF2.ix[sorted(shared_index,reverse=True)]

    merged_df = shared_DF1.join(shared_DF2.loc[:, ['T_mut_per',
 'T_ref_cov',
 'T_mut_cov',
 'N_mut_per',
 'N_ref_cov',
 'N_mut_cov',]],lsuffix='_'+samples[0],rsuffix='_'+samples[1])

    former_cols = ['Chr',
 'Start',
 'End',
 'Ref',
 'Alt',
                   #5
 'Func.refGene',
 'Gene.refGene',
 'GeneDetail.refGene',
 'ExonicFunc.refGene',
 'AAChange.refGene',
 'snp138',
 'CLINSIG',
 'CLNDBN',
 'CLNACC',
 'CLNDSDB',
 'CLNDSDBID']

    processed = ['T_mut_per',
     'T_ref_cov',
     'T_mut_cov',
     'N_mut_per',
     'N_ref_cov',
     'N_mut_cov', ]

    _temp = [_i+'_'+samples[0] for _i in processed] + [_i+'_'+samples[1] for _i in processed]
    _temp.reverse()
    for _i in _temp:
        former_cols.insert(5,_i)
    merged_df.loc[:,former_cols].to_csv(output,index=False)