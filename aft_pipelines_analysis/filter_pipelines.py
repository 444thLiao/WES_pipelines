from genes_list import snp_138_common_init,formatt_col
import pandas,time
import argparse
import os, sys
import pandas as pd
from tqdm import tqdm
from glob import glob
def filter_pipelines2(normal_germline,normal_somatic, tumor_somatic,pair_somatic,output_path,range_bed,pp=[0,2,3,4,5,6],snp_path='/home/liaoth/data/humandb/snp138Common.name.txt'):
    '''
    created at 2017-06-26 advances XK filter pipelines

    :param normal_germline:
    :param tumor_somatic:
    :param pair_somatic:
    :param output_path:
    :param filter_pipeliens: select the filter you want to use,order isn't matter.
    :return:
    '''
    t1 = time.time()
    snp138_common = snp_138_common_init(snp_path)
    TUMOR_S = pd.read_csv(tumor_somatic, index_col=None)
    TUMOR_P = pd.read_csv(pair_somatic, index_col=None)
    NORMAL_S = pd.read_csv(normal_somatic, index_col=None)
    NORMAL_germline = pd.read_csv(normal_germline, index_col=None)
    t2 = time.time()
    print('loading all needed data, using %d s' % (t2-t1))

    TUMOR_S_index = [';'.join([str(_i) for _i in list(TUMOR_S.iloc[_idx, :5])]) for _idx in list(TUMOR_S.index)]
    TUMOR_P_index = [';'.join([str(_i) for _i in list(TUMOR_P.iloc[_idx, :5])]) for _idx in list(TUMOR_P.index)]
    NORMAL_S_index = [';'.join([str(_i) for _i in list(NORMAL_S.iloc[_idx, :5])]) for _idx in list(NORMAL_S.index)]
    NORMAL_germline_index = [';'.join([str(_i) for _i in list(NORMAL_germline.iloc[_idx, :5])]) for _idx in list(NORMAL_germline.index)]

    TUMOR_S.index = TUMOR_S_index
    TUMOR_P.index = TUMOR_P_index
    NORMAL_germline.index = NORMAL_germline_index
    NORMAL_S.index = NORMAL_S_index

    NORMAL_germline_index = pass_filter(NORMAL_germline)
    NORMAL_S_index = F_filter(NORMAL_S,option='higher',threshold=0.03)
    #choose only the germline variants with pass filters and higher than 3% AF.

    #init
    TUMOR_S_filtered_index = TUMOR_S_index
    TUMOR_P_filtered_index = TUMOR_P_index

    descriptions = []
    counts = []
    descriptions.append('Ori')
    counts.append([len(TUMOR_S),len(TUMOR_P)])
    t3 = time.time()
    print('Initing all listed data, using %d s' % (t3-t2))

    if 0 in pp:
        special_S_index = clinvar_filter(TUMOR_S)
        special_P_index = clinvar_filter(TUMOR_P)
        print('finish extracting the clinvar imp var......')
        descriptions.append('Extract actionable ')
        counts.append([len(special_S_index), len(special_P_index)])
    else:
        special_S_index = special_P_index = []

    if 1 in pp:
        TUMOR_S_filtered_index = F_filter(TUMOR_S, option='higher', threshold=0.1)
        TUMOR_P_filtered_index = F_filter(TUMOR_P, option='higher', threshold=0.1)
        print('finish frequency filter.....')
        descriptions.append('Frequency')
        counts.append([len(TUMOR_S_filtered_index),
                       len(TUMOR_P_filtered_index)])

    if 2 in pp:
        TUMOR_S_filtered_index = cov_filter_info_Version(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = cov_filter_info_Version(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print('finish coverage filter.....')
        descriptions.append('Coverage')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 3 in pp:
        TUMOR_S_filtered_index = snp_common(snp138_common, TUMOR_S.loc[TUMOR_S_filtered_index+special_S_index, :],)
        TUMOR_P_filtered_index = snp_common(snp138_common, TUMOR_P.loc[TUMOR_P_filtered_index+special_P_index, :])
        print('finish snp common filter.....')
        descriptions.append('Snp_common')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 4 in pp:
        TUMOR_S_filtered_index = fun_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = fun_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print('finish Gene function filter.....')
        descriptions.append('Gene_Function')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])

    if 5 in pp:
        TUMOR_S_filtered_index = pass_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = pass_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print('finish pass filter.....')
        descriptions.append('Pass_filter')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 6 in pp:
        TUMOR_P_filtered_index = set(TUMOR_P_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        TUMOR_S_filtered_index = set(TUMOR_S_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        print('finish de-germline filter.....')
        descriptions.append('De-germline_filter')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 7 in pp:
        TUMOR_S_remained_index = intarget_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :],range_list = range_bed)
        TUMOR_P_remained_index = intarget_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :],range_list = range_bed)
        print("finish 'target and TAF>=NAF' filter.....")
        descriptions.append('target and TAF>=NAF')
        counts.append([len(TUMOR_S_filtered_index), len(TUMOR_P_filtered_index)])

    result1 = TUMOR_S.loc[sorted(set(TUMOR_S_filtered_index).intersection(set(TUMOR_P_filtered_index))),:]
    # Remove same variants in order to merge.
    result2 = TUMOR_P.loc[sorted(set(TUMOR_P_filtered_index))]
    result = pandas.concat([result1,result2])

    descriptions.append('Deduplication')
    counts.append([len(result)])

    remain_idxs = clinvar_filter(result)
    descriptions.append('clinvar_filter')
    counts.append([len(set(list(result.index)).difference(set(remain_idxs)))])
    with open(output_path,'w') as f1:
        result.loc[set(list(result.index)).difference(set(remain_idxs)),formatt_col].to_csv(f1,index=False)
    tmp = pd.DataFrame(index=descriptions, columns=['Single', 'Paired'])
    with open(output_path.replace('.csv', '.summary'), 'w') as f1:
        tmp.loc[:, 'Single'] = [_[0] for _ in counts]
        tmp.loc[:, 'Paired'] = [_[-1] for _ in counts]
        tmp.to_csv(f1)
    # print(descriptions,counts)



if __name__ == '__main__':
    ############################################################
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', dest='input', type=str,required = True,
                                help="input dir ")
    parse.add_argument('-o', dest='output', type=str,required = True,
                                help="output dir")
    parse.add_argument('-s', dest='setting', type=str,required = True,
                                help="setting file")
    parse.add_argument('-b', dest='bed', type=str,required = True,
                                help="bed file for WES")
    args = parse.parse_args()
    input_dir = os.path.abspath(args.input)
    setting_file = args.setting
    output_dir = os.path.abspath(args.output)
    bed_file = os.path.abspath(args.bed)

    BED_INFO = pd.read_csv(bed_file, sep='\t', header=None,
                           index_col=None)
    range_list = []
    for idx in tqdm(range(BED_INFO.shape[0])):
        s, e = BED_INFO.iloc[idx, [1, 2]].values.tolist()
        range_list += range(s, e + 1)

    ############################################################
    import_path = os.path.dirname(setting_file)
    sys.path.insert(0, import_path)
    from setting import *
    from filters import *
    ############################################################
    samples_ids = [csv_f.split(".merged.anno.")[0] for csv_f in glob(os.path.join(input_dir,'**','*.csv'),
                                                                     recursive=True)]
    pair_name = [i for i in samples_ids if NORMAL_SIG not in i and TUMOR_SIG not in i]
    for each_pair in pair_name:
        if each_pair.count('-') == 2:
            normal_name = each_pair.rpartition('-')[0] + NORMAL_SIG
            tumor_name = each_pair.rpartition('-')[0] + TUMOR_SIG + '-' + each_pair.rpartition('-')[2]
        else:
            normal_name = each_pair + NORMAL_SIG
            tumor_name = each_pair + TUMOR_SIG
        germline = glob.glob(os.path.join(input_dir, 'germline', normal_name + '.merged*.csv'))[0]
        somatic_normal = glob.glob(os.path.join(input_dir, 'somatic', normal_name + '.mt2*.csv'))[0]
        somatic_tumor = glob.glob(os.path.join(input_dir, 'somatic', tumor_name + '.mt2*.csv'))[0]
        somatic_pair = glob.glob(os.path.join(input_dir, 'somatic', each_pair + '.mt2*.csv'))[0]
        output_file = os.path.join(output_dir,'%s_all_except_AF_depth.csv' % each_pair)
        if not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        print("filter_pipelines2(%s,%s,%s,%s,%s,pp=[3,4,5,6])" % (germline,somatic_normal,somatic_tumor,somatic_pair,output_file))
        filter_pipelines2(germline,
                          somatic_normal,
                          somatic_tumor,
                          somatic_pair,
                          output_file,
                          range_bed=range_list,
                          pp=[3,4,5,6],
                          snp_path=snp138_common_file)
        filter_pipelines2(germline,
                          somatic_normal,
                          somatic_tumor,
                          somatic_pair,
                          output_file.replace('except_AF_depth.csv','except_AF_depth_PASS.csv'),
                          pp=[3, 4, 6],
                          range_bed=range_list,
                          snp_path=snp138_common_file)
        # print()


    # filter_pipelines2('/home/liaoth/project/170801_XK/result/germline/XK-27W.merged.anno.csv.hg19_multianno.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27W.mt2.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27T_with_info.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27_with_info.csv',
    #                  output_file,
    #                  pp=arg_list)