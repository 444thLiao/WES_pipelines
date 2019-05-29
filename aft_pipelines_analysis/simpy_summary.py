import pandas as pd
from collections import Counter,defaultdict
import os
def make_summary(input_csv):
    count = 0
    if type(input_csv) == str:
        result = pd.DataFrame(columns=['Gene_fun', 'Exonic_fun', 'count'])
        cache_df = pd.read_csv(input_csv)
        funcs = set(cache_df.loc[:,'Func.refGene'])
        for _func in funcs:
            counter_yet = Counter(list(cache_df[cache_df.loc[:,'Func.refGene'] == _func].loc[:,'ExonicFunc.refGene']))
            for _k in sorted(counter_yet.keys()):
                result.loc[count, :] = [_func, _k, str(counter_yet[_k])]
                count += 1
        return result
    else:
        columns = ['Gene_fun', 'Exonic_fun']
        for _csv in input_csv:
            samples_name = os.path.basename(_csv).split('.merged.anno')[0]
            columns.append(samples_name)
        result = pd.DataFrame(columns=columns)
        for _csv in input_csv:
            samples_name = os.path.basename(_csv).split('.merged.anno')[0]
            cache_df = pd.read_csv(_csv)
            counter_yet = Counter(
                [';;;'.join(_) for _ in cache_df.loc[:, ['Func.refGene', 'ExonicFunc.refGene']].values.tolist()])
            for _k in sorted(counter_yet.keys()):
                result.loc[_k, samples_name] = counter_yet[_k]
                result.loc[_k, 'Gene_fun'] = _k.split(';;;')[0]
                result.loc[_k, 'Exonic_fun'] = _k.split(';;;')[1]
        result.index = range(result.shape[0])
        result.sort_values(['Gene_fun', 'Exonic_fun'], inplace=True)
        result.fillna(0, inplace=True)
        result.loc['total',:] = result.sum(0)
        result.loc['total','Gene_fun'] ='total'
        result.loc['total', 'Exonic_fun'] = ''
        return result


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
    biomarkers = pd.read_csv(snvs_indels_biomarkers, sep='\t')
    tiers = pd.read_csv(snvs_indels_tiers, sep='\t')
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
    tiers = pd.read_csv(snvs_indels_tiers, sep='\t')

    tiers.index = tiers.GENOMIC_CHANGE
    cancer_hotspots = tiers[tiers.TIER == 'TIER 2'].loc[:,cols]
    cancer_hotspots.rename(columns={'SYMBOL': 'GENE'}, inplace=True)
    with open(output_csv,'w') as f1:
        cancer_hotspots.to_csv(f1,index=False)


if __name__ == '__main__':
    # TODO
    # finish aft pipelines
    pass