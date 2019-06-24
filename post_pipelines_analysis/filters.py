import pandas as pd
from tqdm import tqdm

def extract2dict(Otherinfo):
    _info = Otherinfo.split('\t')
    _n = _info[-2].split(':')
    _v = _info[-1].split(':')
    _query = dict(zip(_n, _v))

    extrat_dict = _info[7].split(';')
    extrat_dict = dict([_.split('=') for _ in extrat_dict if '=' in _])
    _query.update(extrat_dict)
    return _query
def get_af(_info):
    " get depth frequency of alt "
    _query = extract2dict(_info)

    af_info = _query.get('AF', '0')
    if ',' in af_info:
        af_info = float(af_info.split(',')[0])
    else:
        af_info = float(af_info)
    depth_cal_af = _query.get('AD',)
    try:
        ref_d,alt_d = map(float,depth_cal_af.split(','))
    except:
        import pdb;pdb.set_trace()
    depth_af = alt_d / (ref_d+alt_d) * 100

    return depth_af

def snp_common(snpfile_list,file_df):
    "return not common SNP"
    common_snp_index = set(file_df[file_df.snp138.isin(snpfile_list)].index)
    rare_snp = list(set(file_df.index).difference(common_snp_index))

    return rare_snp

def pass_filter(file_df,match=["PASS"] ):
    if 'Otherinfo' not in file_df.columns:
        return 'Wrong df'
    if type(match) != set:
        match = set(match)
    pass_bucket = []
    for idx,row in file_df.iterrows():
        _info = row['Otherinfo']
        _query = str(_info).split('\t')
        if len(_query) <= 6:
            # continue
            raise Exception("too short")
        try:

            if match.intersection(set(_query[6].split(';'))):
                pass_bucket.append(idx)
        except:
            print(idx,'it should be wrong')
    return pass_bucket



def cov_filter_info_Version(file_df,threshold = 10):
    if 'N_mut_cov' not in file_df:
        return 'Wrong Columns'
    file_df = file_df.loc[file_df.N_mut_cov.astype(str) != 'Off target', :]
    file_df = file_df.loc[file_df.T_mut_cov.astype(str) != 'Off target', :]
    _N_cov = file_df.loc[:, ['N_ref_cov', 'N_mut_cov']].astype(float).sum(1)
    _T_cov = file_df.loc[:, ['T_ref_cov', 'T_mut_cov']].astype(float).sum(1)
    high_cov_bucket = list(file_df.loc[(_N_cov >= threshold) & (_T_cov >= threshold) & (
                file_df.loc[:, 'T_mut_per'] > file_df.loc[:, 'N_mut_per']), :].index)

    return high_cov_bucket


def F_filter(file_df,option='lower',threshold=0.03):
    if 'Otherinfo' not in file_df.columns:
        return 'Wrong df'

    F_bucket = []
    for idx,row in file_df.iterrows():
        _info = row['Otherinfo']
        _query = extract2dict(_info)

        af_info = _query.get('AF','0')
        if ',' in af_info:
            af_info = float(af_info.split(',')[0])
        else:
            af_info = float(af_info)
        if option == 'higher':
            if af_info >= threshold:
                F_bucket.append(idx)
        elif option == 'lower':
            if af_info <= threshold:
                F_bucket.append(idx)
    return F_bucket


def intarget_filter(file_df,range_list):
    "return SNP in sequencing range"
    subset_df = file_df.loc[file_df.loc[:, 'Start'].isin(range_list), :]
    return list(subset_df.index)

def clinvar_filter(file_df,not_in = ('benign','likely benign')):
    clin_imp_bucket = []
    for _index in list(file_df.index):
        _info = str(file_df.loc[_index,'CLINSIG']).lower()
        _info = _info.split('|')
        if _info != ['.'] and not set(not_in).intersection(set(_info)):
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
                print('get a genes values = '+str(file_df.loc[_idx,'Gene.refGene'])+" which can't be split")
                continue
            for _g in genes:
                if _g in _in and _idx not in gene_var_bucket:
                    gene_var_bucket.append(_idx)
    return gene_var_bucket

def gene_filter2(file_path,coverage_info):
    pieces_index = []
    file_df = pd.read_csv(file_path,index_col=False)
    coverage_df = pd.read_csv(coverage_info,index_col=0)
    for _chr in list(set(file_df.loc[:,'Chr'])):
        pieces_index += list(file_df[file_df.loc[:,'Start'].isin(list(coverage_df[coverage_df.loc[:,'Chr']==_chr].loc[:,'Posistion']))].index)
    return list(set(pieces_index))