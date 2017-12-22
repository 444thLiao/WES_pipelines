from pandas import DataFrame as df
from genes_list import snp_138_common_init,formatt_col
from filters import *
import pandas,time
import argparse


def filter_pipelines(normal_germline,normal_somatic, tumor_somatic,pair_somatic,output_path,pp=[0,1,2,3,4,5,6]):
    '''
    created at 2017-05-09 advances XK filter pipelines

    :param normal_germline:
    :param tumor_somatic:
    :param pair_somatic:
    :param output_path:
    :param filter_pipeliens: select the filter you want to use,order isn't matter.
    :return:
    '''
    t1 = time.time()
    snp138_common = snp_138_common_init()
    TUMOR_S = df.from_csv(tumor_somatic,index_col=False)
    TUMOR_P = df.from_csv(pair_somatic,index_col=False)
    NORMAL_S = df.from_csv(normal_somatic,index_col=False)
    NORMAL_germline = df.from_csv(normal_germline,index_col=False)
    t2 = time.time()
    print 'loading all needed data, using %d s' % (t2-t1)
    TUMOR_S_index = [';'.join(_i.split('\t')[:5]) for _i in list(TUMOR_S.Otherinfo)]
    TUMOR_P_index = [';'.join(_i.split('\t')[:5]) for _i in list(TUMOR_P.Otherinfo)]
    NORMAL_S_index = [';'.join(_i.split('\t')[:5]) for _i in list(NORMAL_S.Otherinfo)]
    NORMAL_germline_index = [';'.join(_i.split('\t')[:5]) for _i in list(NORMAL_germline.Otherinfo)]

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

    descritions = []
    counts = []
    t3 = time.time()
    print 'Initing all listed data, using %d s' % (t3-t2)

    if 0 in pp:
        special_S_index = clinvar_filter(TUMOR_S)
        special_P_index = clinvar_filter(TUMOR_P)
        print 'finish extracting the clinvar imp var......'
    else:
        special_S_index = special_P_index = []


    if 1 in pp:
        TUMOR_S_filtered_index = F_filter(TUMOR_S, option='higher', threshold=0.1)
        TUMOR_P_filtered_index = F_filter(TUMOR_P, option='higher', threshold=0.1)
        print 'finish frequency filter.....'
    descritions.append('Frequency')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 2 in pp:
        TUMOR_S_filtered_index = cov_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = cov_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish coverage filter.....'
    descritions.append('Coverage')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 3 in pp:
        TUMOR_S_filtered_index = snp_common(snp138_common, TUMOR_S.loc[TUMOR_S_filtered_index+special_S_index, :])
        TUMOR_P_filtered_index = snp_common(snp138_common, TUMOR_P.loc[TUMOR_P_filtered_index+special_P_index, :])
        print 'finish snp common filter.....'
    descritions.append('Snp_common')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 4 in pp:
        TUMOR_S_filtered_index = fun_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = fun_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish Gene function filter.....'
    descritions.append('Gene_Function')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 5 in pp:
        TUMOR_S_filtered_index = pass_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = pass_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish pass filter.....'
    descritions.append('Pass_filter')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 6 in pp:
        TUMOR_P_filtered_index = set(TUMOR_P_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        TUMOR_S_filtered_index = set(TUMOR_S_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        print 'finish de-germline filter.....'
    descritions.append('De-germline_filter')
    counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])


    result1 = TUMOR_S.loc[sorted(set(TUMOR_S_filtered_index).difference(set(TUMOR_P_filtered_index))),:]
    # Remove same variants in order to merge.
    result2 = TUMOR_P.loc[sorted(set(TUMOR_P_filtered_index))]
    result = pandas.concat([result1,result2])

    descritions.append('Deduplication')
    counts.append([len(result)])

    remain_idxs = clinvar_filter(result)

    descritions.append('clinvar_filter')
    counts.append([len(set(list(result.index)).difference(set(remain_idxs)))])
    with open(output_path,'w') as f1:
        result.loc[set(list(result.index)).difference(set(remain_idxs)),:].to_csv(f1,index=False)
    print descritions,counts


def filter_pipelines2(normal_germline,normal_somatic, tumor_somatic,pair_somatic,output_path,pp=[0,2,3,4,5,6]):
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
    snp138_common = snp_138_common_init()
    TUMOR_S = df.from_csv(tumor_somatic,index_col=False)
    TUMOR_P = df.from_csv(pair_somatic,index_col=False)
    NORMAL_S = df.from_csv(normal_somatic,index_col=False)
    NORMAL_germline = df.from_csv(normal_germline,index_col=False)
    t2 = time.time()
    print 'loading all needed data, using %d s' % (t2-t1)

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

    descritions = []
    counts = []
    descritions.append('Ori')
    counts.append([len(TUMOR_S),len(TUMOR_P)])
    t3 = time.time()
    print 'Initing all listed data, using %d s' % (t3-t2)

    if 0 in pp:
        special_S_index = clinvar_filter(TUMOR_S)
        special_P_index = clinvar_filter(TUMOR_P)
        print 'finish extracting the clinvar imp var......'
        descritions.append('Extract actionable ')
        counts.append([len(special_S_index), len(special_P_index)])
    else:
        special_S_index = special_P_index = []

    if 1 in pp:
        TUMOR_S_filtered_index = F_filter(TUMOR_S, option='higher', threshold=0.1)
        TUMOR_P_filtered_index = F_filter(TUMOR_P, option='higher', threshold=0.1)
        print 'finish frequency filter.....'
        descritions.append('Frequency')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])

    if 2 in pp:
        TUMOR_S_filtered_index = cov_filter_info_Version(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = cov_filter_info_Version(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish coverage filter.....'
        descritions.append('Coverage')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 3 in pp:
        TUMOR_S_filtered_index = snp_common(snp138_common, TUMOR_S.loc[TUMOR_S_filtered_index+special_S_index, :])
        TUMOR_P_filtered_index = snp_common(snp138_common, TUMOR_P.loc[TUMOR_P_filtered_index+special_P_index, :])
        print 'finish snp common filter.....'
        descritions.append('Snp_common')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 4 in pp:
        TUMOR_S_filtered_index = fun_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = fun_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish Gene function filter.....'
        descritions.append('Gene_Function')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])

    if 5 in pp:
        TUMOR_S_filtered_index = pass_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = pass_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print 'finish pass filter.....'
        descritions.append('Pass_filter')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])
    if 6 in pp:
        TUMOR_P_filtered_index = set(TUMOR_P_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        TUMOR_S_filtered_index = set(TUMOR_S_filtered_index).difference(set(NORMAL_germline_index),set(NORMAL_S_index))
        print 'finish de-germline filter.....'
        descritions.append('De-germline_filter')
        counts.append([len(TUMOR_S_filtered_index),len(TUMOR_P_filtered_index)])


    result1 = TUMOR_S.loc[sorted(set(TUMOR_S_filtered_index).difference(set(TUMOR_P_filtered_index))),:]
    # Remove same variants in order to merge.
    result2 = TUMOR_P.loc[sorted(set(TUMOR_P_filtered_index))]
    result = pandas.concat([result1,result2])

    descritions.append('Deduplication')
    counts.append([len(result)])

    remain_idxs = clinvar_filter(result)
    #import pdb;pdb.set_trace()
    descritions.append('clinvar_filter')
    counts.append([len(set(list(result.index)).difference(set(remain_idxs)))])
    with open(output_path,'w') as f1:
        result.loc[set(list(result.index)).difference(set(remain_idxs)),formatt_col].to_csv(f1,index=False)
    print descritions,counts


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-pp',dest='pipelines', type=str,
                    help="Pipelines number specify by ','. ")
    parse.add_argument('-o', dest='output', type=str,required = True,
                                help="output file ")
    args = parse.parse_args()
    arg_list = [int(_i) for _i in args.pipelines.split(',')]
    output_file = args.output
    # filter_pipelines2('/home/liaoth/project/170602_XK/server_result/germline/XK-8W_S18.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-8-2/XK-8W.mt2.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-8-2/XK-8T-2_with_info.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-8-2/XK-8-2_with_info.csv',
    #                  output_file,
    #                  pp=arg_list)
    # filter_pipelines2('/home/liaoth/project/170602_XK/server_result/germline/XK-8W_S18.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-8W.mt2.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/result/XK-8/FOR_COMPARE/170123/XK-8T.mt2_with_infos.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/result/XK-8/FOR_COMPARE/170123/XK-8.mt2_with_infos.csv',
    #                  output_file,
    #                  pp=arg_list)
    # filter_pipelines2('/home/liaoth/project/170602_XK/server_result/germline/XK-2W.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-2/XK-2W.mt2.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-2/XK-2T-2_with_info.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-2/XK-2-2_with_info.csv',
    #                  output_file,
    #                  pp=arg_list)

    # filter_pipelines2('/home/liaoth/project/170602_XK/server_result/germline/XK-2W.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/somatic/XK-2W.mt2.merged.anno.hg19_multianno.csv',
    #                   '/home/liaoth/project/170602_XK/server_result/result/XK-2/FOR_COMPARE/XK-2T_S20.mt2_with_infos.csv',
    #                  '/home/liaoth/project/170602_XK/server_result/result/XK-2/FOR_COMPARE/XK-2.mt2_with_infos.csv',
    #                  output_file,
    #                  pp=arg_list)
    filter_pipelines2('/home/liaoth/project/170801_XK/result/germline/XK-27W.merged.anno.csv.hg19_multianno.csv',
                     '/home/liaoth/project/170801_XK/result/somatic/XK-27W.mt2.merged.anno.hg19_multianno.csv',
                     '/home/liaoth/project/170801_XK/result/somatic/XK-27T_with_info.csv',
                     '/home/liaoth/project/170801_XK/result/somatic/XK-27_with_info.csv',
                     output_file,
                     pp=arg_list)
