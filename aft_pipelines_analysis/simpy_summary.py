from pandas import DataFrame as df
from collections import Counter
def make_summary(input_csv):
    result = 'Gene_fun\tExonic_fun\tcount\n'
    cache_df = df.from_csv(input_csv,index_col=False)
    funcs = set(cache_df.loc[:,'Func.refGene'])
    for _func in funcs:
        counter_yet = Counter(list(cache_df[cache_df.loc[:,'Func.refGene'] == _func].loc[:,'ExonicFunc.refGene']))
        for _k in sorted(counter_yet.keys()):
            result += '%s\t%s\t%s\n' % (_func,_k,str(counter_yet[_k]))
    return result



def extract_unpass_snv(withpass_filter_csv,withoutpass_filter_csv):
    withpass_filter_df = df.from_csv(withpass_filter_csv, index_col=False)
    withoutpass_filter_df = df.from_csv(withoutpass_filter_csv, index_col=False)

    withpass_filter_index = [';'.join([str(_i) for _i in list(withpass_filter_df.iloc[_idx, :5])]) for _idx in list(withpass_filter_df.index)]
    withoutpass_filter_index = [';'.join([str(_i) for _i in list(withoutpass_filter_df.iloc[_idx, :5])]) for _idx in list(withoutpass_filter_df.index)]

    withpass_filter_df.index = withpass_filter_index
    withoutpass_filter_df.index = withoutpass_filter_index

    result_index = set(withoutpass_filter_index).difference(set(withpass_filter_index))
    return withoutpass_filter_df.loc[result_index,:]



def construct_biomarkers_csv(snvs_indels_biomarkers,snvs_indels_tiers,output_csv):
    cols = ['GENOMIC_CHANGE',
            'GENOME_VERSION',
            'PROTEIN_CHANGE',
            'PROTEIN_DOMAIN',
            'CDS_CHANGE',
            'SYMBOL',
            'CONSEQUENCE',
            'BM_CLINICAL_SIGNIFICANCE',
            'BM_EVIDENCE_LEVEL',
            'BM_EVIDENCE_TYPE',
            'BM_EVIDENCE_DIRECTION',
            'BM_DISEASE_NAME',
            'BM_DRUG_NAMES',
            'BM_CITATION',]
    biomarkers = df.from_csv(snvs_indels_biomarkers,sep='\t',index_col=False)
    tiers = df.from_csv(snvs_indels_tiers,sep='\t',index_col=False)
    biomarkers.index = biomarkers.GENOMIC_CHANGE
    tiers.index = tiers.GENOMIC_CHANGE

    new_biomarkers = biomarkers.join(tiers.loc[:,['PROTEIN_CHANGE','PROTEIN_DOMAIN','CDS_CHANGE']]).loc[:,cols]
    new_biomarkers.rename(columns = {'SYMBOL':'GENE'},inplace=True)
    with open(output_csv,'w') as f1:
        new_biomarkers.to_csv(f1,index=False)




def extract_tier2(snvs_indels_tiers,output_csv):
    cols = ['GENOMIC_CHANGE',
            'GENOME_VERSION',
            'DBSNP',
            'PROTEIN_CHANGE',
            'PROTEIN_DOMAIN',
            'CDS_CHANGE',
            'SYMBOL',
            'CONSEQUENCE',
            'BM_CLINICAL_SIGNIFICANCE',
            'BM_EVIDENCE_LEVEL',
            'BM_EVIDENCE_TYPE',
            'BM_EVIDENCE_DIRECTION',
            'BM_DISEASE_NAME',
            'BM_DRUG_NAMES',
            'BM_CITATION', ]
    tiers = df.from_csv(snvs_indels_tiers,sep='\t',index_col=False)

    tiers.index = tiers.GENOMIC_CHANGE
    cancer_hotspots = tiers[tiers.TIER == 'TIER 2'].loc[:,cols]
    cancer_hotspots.rename(columns={'SYMBOL': 'GENE'}, inplace=True)
    with open(output_csv,'w') as f1:
        cancer_hotspots.to_csv(f1,index=False)