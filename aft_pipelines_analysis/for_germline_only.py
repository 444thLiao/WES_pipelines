import os
from glob import glob
import filters
import pandas as pd
from tqdm import tqdm
from genes_list import snp_138_common_init

odir = "/home/liaoth/data2/project/ZHJ_WES/filtered_o"
indir = "/home/liaoth/data2/project/ZHJ_WES/raw_o"
bed_file = "/home/liaoth/data2/project/ZHJ_WES/Sureselect_V6_COSMIC_formal.bed"
############################################################
snp_path = '/home/liaoth/data/humandb/snp138Common.name.txt'
normal_germline = '/home/liaoth/data2/project/ZHJ_WES/output/HUA.merged.anno.hg19_multianno.csv'
snp138_common = snp_138_common_init(snp_path)

BED_INFO = pd.read_csv(bed_file, sep='\t', header=None,
                       index_col=None)
range_list = []
# for idx in tqdm(range(BED_INFO.shape[0])):
for idx, row in tqdm(BED_INFO.iterrows(),
                     total=BED_INFO.shape[0]):
    s, e = row[1], row[2]
    #
    # s, e = BED_INFO.iloc[idx, [1, 2]].values.tolist()
    range_list += range(s, e + 1)

remained_cols = ['Func.refGene',
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
for normal_germline in tqdm(glob(os.path.join(indir, '*.merged.anno.*.csv'))):

    NORMAL_germline = pd.read_csv(normal_germline, index_col=None)
    NORMAL_germline_index = [';'.join([str(_i) for _i in list(NORMAL_germline.iloc[_idx, :5])]) for _idx in list(NORMAL_germline.index)]

    NORMAL_germline.loc[:,'AF'] = list(map(filters.get_af,list(NORMAL_germline.Otherinfo)))
    pass_index = filters.pass_filter(NORMAL_germline, match={"PASS"})
    intarget_index = filters.intarget_filter(NORMAL_germline, range_list=range_list)
    basic_snp = set.intersection(set(pass_index),
                                 set(intarget_index))

    clinvar_index = filters.clinvar_filter(NORMAL_germline.loc[basic_snp, :])

    high_freq_snp = filters.F_filter(NORMAL_germline.loc[basic_snp, :], option='higher', threshold=0.1)

    rare_SNP_index = filters.snp_common(snp138_common, NORMAL_germline.loc[basic_snp, :], )
    fun_index = filters.fun_filter(NORMAL_germline.loc[basic_snp, :])

    ############################################################
    sample_id = os.path.basename(normal_germline).split('.merged.anno')[0]
    for filtered_t, idx in zip(["clinvar", "high_freq", "rare_SNP", "functional"],
                               [clinvar_index, high_freq_snp, rare_SNP_index, fun_index]):
        os.makedirs(os.path.join(odir, sample_id), exist_ok=True)
        base_name = os.path.join(odir, sample_id, "%s_%s_filtered.csv" % (sample_id, filtered_t))

        filtered_df = NORMAL_germline.loc[idx, remained_cols]
