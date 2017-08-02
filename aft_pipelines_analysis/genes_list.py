from pandas import DataFrame as df

def snp_138_common_init(path = '/home/liaoth/data/humandb/snp138Common.name.txt'):
    '''

    :param path: snp 138 common file
    :return: list of every common snp.
    '''
    snp138_common = open(path).readline()
    snp138_common = [_i.replace('\n', '') for _i in snp138_common]
    return snp138_common

def ClearSeq_pos_init(path = '/home/liaoth/project/170123_NY/ClearSeq.bed'):
    '''

    :param path: ClearSeq 151 gene info file
    :return: pandas DataFrame with Start,End info. Using Utils.construct_pos_list() to create pos dict.
    '''
    gene_151 = df.from_csv(path, sep='\t', index_col=False)
    return gene_151



sum_sig_genes = ['AKT1', 'AKT2', 'AKT3', 'ARID1A', 'ARID1B', 'ARID2', 'ASCL4', 'ATM', 'BRAF', 'CDKN2A', 'COBL',
                 'CREBBP', 'CTNNB1', 'CUL3', 'EGFR', 'EP300', 'EPHA7', 'ERBB2', 'ERBB3', 'FGFR1', 'FGFR2', 'FGFR3',
                 'FOXP2', 'HRAS', 'KEAP1', 'KMT2D', 'KRAS', 'MAP2K1', 'MET', 'MGA', 'MLL', 'NF1', 'NFE2L2', 'NOTCH1',
                 'NOTCH2', 'NRAS', 'PIK3CA', 'PTEN', 'RASA1', 'RB1', 'RBM10', 'RIT1', 'SETD2', 'SLIT2', 'SMAD4',
                 'SMARCA4', 'SOX2', 'STK11', 'TP53', 'TP63', 'TSC1', 'TSC2', 'U2AF1','ALK','DDR2']
##Lung Cancer relative important genes.



formatt_col = ['Chr',
 'Start',
 'End',
 'Ref',
 'Alt',
 'T_mut_per',
 'T_ref_cov',
 'T_mut_cov',
 'N_mut_per',
 'N_ref_cov',
 'N_mut_cov',
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