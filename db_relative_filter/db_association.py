
from pandas import DataFrame as df
import re


def filter_base_pos(input_csv,
                    input_query,
                    output,
                    q_csv = '\t',
                    q_cols = ['chromosome','start','stop'],
                    cols_prefix='test'):
    subject_csv = df.from_csv(input_csv,index_col=False)
    query_csv = df.from_csv(input_query,sep=q_csv,index_col=False)

    exist_pair = []

    for _idx in list(query_csv.index):
        q_chr = query_csv.loc[_idx,q_cols[0]]
        q_start = query_csv.loc[_idx,q_cols[1]]
        q_end = query_csv.loc[_idx,q_cols[2]]

        s_cache = subject_csv[subject_csv.loc[:,'Chr'] == 'chr'+str(q_chr)]

        start_result = s_cache[s_cache.loc[:,'Start'] == q_start]

        if len(start_result) > 0:
            exist_pair.append((list(start_result[start_result.End == q_end].index),_idx))
    if exist_pair:
        with open(output,'w') as f1:
            try:
                s_idxs = [_i[0][0] for _i in exist_pair if len(_i[0])==1]
                q_idxs = [_i[1] for _i in exist_pair if len(_i[0])==1]
                subject_csv.loc[s_idxs, '%s_db' % cols_prefix] = [';'.join([str(_k) for _k in list(query_csv.loc[_j,:])]) for _j in q_idxs]
            except:
                import pdb;pdb.set_trace()

            subject_csv.loc[s_idxs,:].to_csv(f1,index=False)
    return exist_pair

def filter_base_alteration(input_csv,
                           input_query,
                           output,
                           q_col = 'Alteration',
                           s_col = 'AAChange.refGene',
                           cols_prefix = 'test',
                           altera_pattern = '[A-Z]\d+[A-Z]?'):
    subject_csv = df.from_csv(input_csv,index_col=False)
    query_csv = df.from_csv(input_query,sep='\t',index_col=False)

    subject_csv = subject_csv[subject_csv.loc[:,s_col] != '.']

    exist_pair_bucket = []

    for _idx in list(query_csv.index):
        cache = subject_csv[subject_csv.loc[:,'Gene.refGene'] == query_csv.loc[_idx,'Gene']]
        for _sidx in list(cache.index):
            if query_csv.loc[_idx,q_col].replace('*','') in cache.loc[_sidx,s_col] and re.findall(altera_pattern,query_csv.loc[_idx,q_col].replace('*','')):
                # remove * in some of the alteration.

                #print query_csv.loc[_idx,:],cache.loc[_sidx,:]
                exist_pair_bucket.append((_sidx,_idx)) # (s,q)

    with open(output,'w') as f1:
        cache = subject_csv.loc[[_s for _s,_q in exist_pair_bucket],:]
        for _s, _q in exist_pair_bucket:
            cache.loc[_s, '%s_db' % cols_prefix] = ';'.join([str(_j) for _j in list(query_csv.loc[_q, :])])
        cache.to_csv(f1,index=False)
    return exist_pair_bucket



def filter_base_alteration2(input_csv,
                            input_query,
                            output,
                            q_col = 'Alteration',
                            s_col = 'AAChange.refGene',cols_prefix = 'test',altera_pattern = '[A-Z]+\d+[A-Z]+'):
    #for  '/home/liaoth/data/Cancer_db/cancer_genome_interpreter/cgi_biomarkers_20170208/cgi_biomarkers_per_variant.tsv'
    subject_csv = df.from_csv(input_csv,index_col=False)
    query_csv = df.from_csv(input_query,sep='\t',index_col=False)

    subject_csv = subject_csv[subject_csv.loc[:,s_col] != '.']

    exist_pair_bucket = []

    for _idx in list(query_csv.index):
        gene = query_csv.loc[_idx,q_col].split(':')[0] if ':' in str(query_csv.loc[_idx,q_col]) else query_csv.loc[_idx,q_col].split('__')[0]
        cache = subject_csv[subject_csv.loc[:,'Gene.refGene'] == gene]
        for _sidx in list(cache.index):
            if re.findall(altera_pattern, query_csv.loc[_idx, q_col].replace('*', 'X')):
                for each in re.findall(altera_pattern, query_csv.loc[_idx, q_col].replace('*', 'X')):
                    if each in cache.loc[_sidx,s_col] and each not in gene:
                        exist_pair_bucket.append((_sidx, _idx))
                        print each,cache.loc[_sidx,s_col],gene
                        break
    with open(output,'w') as f1:
        cache = subject_csv.loc[[_s for _s,_q in exist_pair_bucket],:]
        for _s, _q in exist_pair_bucket:
            cache.loc[_s, '%s_db' % cols_prefix] = ';'.join([str(_j) for _j in list(query_csv.loc[_q, :])])
        cache.to_csv(f1,index=False)
    return exist_pair_bucket


with open('/home/liaoth/project/NY_project_shanghai/from_server/combinded_hg19_cosmic_v77.vcf') as f1:
    cosmic = [_.replace('\n', '') for _ in f1.readlines() if 'chr' in _]
    cosmic_pos = set([';'.join(_.split('\t')[:2]) for _ in cosmic])

def filter_out_cosmic(csv_path):
    _cache = df.from_csv(csv_path,index_col=False)
    _cache_list = [';'.join([str(_) for _ in list(_cache.iloc[_i,:2])]) for _i in list(_cache.index)]
    _cache.index = _cache_list
    bucket = []
    for _ in _cache_list:
        if _ in cosmic_pos:
            bucket.append(_)
    return _cache.loc[bucket,:]