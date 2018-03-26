from Utils import extract2dict
from pandas import DataFrame as df
import tqdm

BED_INFO = df.from_csv('/home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed', sep='\t', header=None,
                       index_col=None)
range_list = []
for idx in tqdm.tqdm(range(BED_INFO.shape[0])):
    s, e = BED_INFO.iloc[idx, [1, 2]].values.tolist()
    range_list += range(s, e + 1)

def snp_common(snpfile_list,file_df):

    common_snp_index = set(file_df[file_df.snp138.isin(snpfile_list)].index)
    rare_snp = list(set(file_df.index).difference(common_snp_index))

    return rare_snp

def pass_filter(file_df,_in = ['PASS']):
    if 'Otherinfo' not in file_df:
        return 'Wrong Columns'

    pass_bucket = []
    for _index in list(file_df.index):
        _info = file_df.loc[_index, 'Otherinfo']
        if type(_info) != str:
            import pdb;pdb.set_trace()
        _query = _info.split('\t')[6]
        try:
            if set(_in).intersection(set(_query.split(';'))):
                pass_bucket.append(_index)
        except:
            print _index,'it should be wrong'
    return pass_bucket



def cov_filter_info_Version(file_df,threshold = 10):
    if 'N_mut_cov' not in file_df:
        return 'Wrong Columns'
    file_df.index = range(file_df.shape[0])
    file_df = file_df.loc[file_df.N_mut_cov != 'Off target', :]
    _N_cov = file_df.loc[:, ['N_ref_cov', 'N_mut_cov']].astype(float).sum(1)
    _T_cov = file_df.loc[:, ['T_ref_cov', 'T_mut_cov']].astype(float).sum(1)
    high_cov_bucket = list(file_df.loc[(_N_cov >= threshold) & (_T_cov >= threshold) & (
                file_df.loc[:, 'T_mut_per'] > file_df.loc[:, 'N_mut_per']), :].index)

    return high_cov_bucket


def F_filter(file_df,option='lower',threshold=0.03):
    if 'Otherinfo' not in file_df:
        return 'Wrong Columns'

    lowF_bucket = []
    for _index in list(file_df.index):
        _info = file_df.loc[_index, 'Otherinfo']
        _query = extract2dict(_info)

        af_info = _query['AF']
        if ',' in af_info:
            af_info = float(af_info.split(',')[0])
        else:
            af_info = float(af_info)
        if option == 'higher':
            if af_info >= threshold:
                lowF_bucket.append(_index)
        elif option == 'lower':
            if af_info <= threshold:
                lowF_bucket.append(_index)
    return lowF_bucket


def offtarget_filter(file_df):
    subset_df = file_df.loc[~file_df.loc[:, 'Start'].isin(range_list), :]
    return list(subset_df.index)

def clinvar_filter(file_df,not_in = ['benign','likely benign']):
    clin_imp_bucket = []
    for _index in list(file_df.index):
        _info = file_df.loc[_index,'CLINSIG'].lower()
        _info = _info.split('|')
        if _info != ['.'] and set(not_in).intersection(set(_info)):
            clin_imp_bucket.append(_index)
    return clin_imp_bucket


def fun_filter(file_df,keep = {'Func.refGene':['exonic','exonic;splicing','splicing']},drop = {'ExonicFunc.refGene':['synonymous SNV']}):
    '''

    :param file_df:
    :param keep: dict key is col_name,value is a list of pos_value need to keep.
    :param drop: dict key is col_name,value is a list of pos_value need to drop.
    :return: index list after filtration .
    '''
    if keep:
        for i in keep.keys():
            filter_file = file_df[file_df.loc[:, i].isin(keep[i])]
    else:
        filter_file = file_df
    ####drop
    if drop:
        for i in drop.keys():
            for k in drop[i]:
                filter_file = filter_file[filter_file.loc[:, i] != k]
    return list(filter_file.index)


def filter_according_gene(file_df):
    from .genes_list import sum_sig_genes,ClearSeq_pos_init
    from .Utils import construct_pos_list
    indexs = file_df.index
    gene_lists = list(file_df.loc[:, 'Gene.refGene'])
    query_for = zip(indexs, gene_lists)

    sig_gene_index = []
    for _i, _v in query_for:
        judge = False
        for _each in _v.split(';'):
            if _each in sum_sig_genes:
                judge = True
        if judge:
            sig_gene_index.append(_i)
    sig_gene_index2 = []
    gene_151_pos = construct_pos_list(ClearSeq_pos_init())
    for _chr in gene_151_pos:
        _cache = file_df[file_df.loc[:,'Chr'] == _chr]
        sig_gene_index2 += list(_cache[_cache.loc[:,'Start'].isin(gene_151_pos[_chr])].index)
    return file_df.loc[set(sig_gene_index + sig_gene_index2), :]


def gene_filter(file_df,_in = []):
    gene_var_bucket = []
    indexs = file_df.index
    for _idx in list(indexs):
        if file_df.loc[_idx,'Gene.refGene']:
            try:
                genes = file_df.loc[_idx,'Gene.refGene'].split(';')
            except:
                print 'get a genes values = '+str(file_df.loc[_idx,'Gene.refGene'])+" which can't be split"
                continue
            for _g in genes:
                if _g in _in and _idx not in gene_var_bucket:
                    gene_var_bucket.append(_idx)
    return gene_var_bucket

def gene_filter2(file_path,coverage_info):
    pieces_index = []
    file_df = df.from_csv(file_path,index_col=False)
    coverage_df = df.from_csv(coverage_info)
    for _chr in list(set(file_df.loc[:,'Chr'])):
        pieces_index += list(file_df[file_df.loc[:,'Start'].isin(list(coverage_df[coverage_df.loc[:,'Chr']==_chr].loc[:,'Posistion']))].index)
    return list(set(pieces_index))