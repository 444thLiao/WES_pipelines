from pandas import DataFrame as df
from collections import Counter,defaultdict
import os
def make_summary(input_csv):
    if type(input_csv) == str:
        result = 'Gene_fun\tExonic_fun\tcount\n'
        cache_df = df.from_csv(input_csv,index_col=False)
        funcs = set(cache_df.loc[:,'Func.refGene'])
        for _func in funcs:
            counter_yet = Counter(list(cache_df[cache_df.loc[:,'Func.refGene'] == _func].loc[:,'ExonicFunc.refGene']))
            for _k in sorted(counter_yet.keys()):
                result += '%s\t%s\t%s\n' % (_func,_k,str(counter_yet[_k]))
        return result
    else:
        result = 'Gene_fun\tExonic_fun'
        fields_bucket = defaultdict(list)
        for path in input_csv:
            cache_df = df.from_csv(path,index_col=False)
            funcs = set(cache_df.loc[:, 'Func.refGene'])
            for _func in funcs:
                fields_bucket[_func] += list(set(cache_df[cache_df.loc[:,'Func.refGene'] == _func].loc[:,'ExonicFunc.refGene']))
            result += '\t' + os.path.basename(path).split('_')[0]
        result+='\n'
        for _k in fields_bucket.keys():
            fields_bucket[_k] = list(set(fields_bucket[_k]))

        for _func in fields_bucket.keys():
            for exonic in fields_bucket[_func]:
                for path in input_csv:
                    cache_df = df.from_csv(path, index_col=False)
                    count = len(cache_df[(cache_df.loc[:,'Func.refGene']== _func) & (cache_df.loc[:,'ExonicFunc.refGene']== exonic)])
                    if not result.endswith('\n'):
                        result += '\t%s' % str(count)
                    else:
                        result += '%s\t%s\t%s' % (_func, exonic, str(count))
                result += '\n'
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