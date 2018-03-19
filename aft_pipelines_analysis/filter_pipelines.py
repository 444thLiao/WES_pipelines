from pandas import DataFrame as df
from genes_list import snp_138_common_init,formatt_col
from filters import *
import pandas,time
import argparse
import os, sys, glob
from Utils import pfn
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
    TUMOR_S = df.from_csv(tumor_somatic, index_col=None)
    TUMOR_P = df.from_csv(pair_somatic, index_col=None)
    NORMAL_S = df.from_csv(normal_somatic, index_col=None)
    NORMAL_germline = df.from_csv(normal_germline, index_col=None)
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
    if 7 in pp:
        TUMOR_S_filtered_index = offtarget_filter(TUMOR_S.loc[TUMOR_S_filtered_index, :])
        TUMOR_P_filtered_index = offtarget_filter(TUMOR_P.loc[TUMOR_P_filtered_index, :])
        print "finish 'target and TAF>=NAF' filter....."
        descritions.append('target and TAF>=NAF')
        counts.append([len(TUMOR_S_filtered_index), len(TUMOR_P_filtered_index)])

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
    import pandas as pd
    tmp = pd.DataFrame(index=descritions, columns=['Single', 'Paired'])
    with open(output_path.replace('.csv', '.summary'), 'w') as f1:
        tmp.loc[:, 'Single'] = [_[0] for _ in counts]
        tmp.loc[:, 'Paired'] = [_[-1] for _ in counts]
        tmp.to_csv(f1)
    print descritions,counts



if __name__ == '__main__':
    ############################################################
    csv_output_dir = '/home/liaoth/data2/project/XK_WES/180309_all/output/'
    setting_file = '/home/liaoth/data2/project/XK_WES/180309_all/setting.py'
    ############################################################
    # parse = argparse.ArgumentParser()
    # parse.add_argument('-o', dest='output', type=str,required = True,
    #                             help="output file ")
    # args = parse.parse_args()
    # output_file = args.output
    ############################################################
    import_path = os.path.dirname(setting_file)
    import_name = os.path.basename(setting_file).replace('.py', '')
    sys.path.insert(0, import_path)
    exec "from %s import *" % import_name

    ############################################################
    import pdb;pdb.set_trace()
    samples_ids = [os.path.basename(i).split('.mt2.merged.')[0] for i in glob.glob(os.path.join(csv_output_dir, 'somatic', '*.csv'))]
    pair_name = [i for i in samples_ids if NORMAL_SIG not in i and TUMOR_SIG not in i]
    for each_pair in pair_name:
        if each_pair.count('-') == 2:
            normal_name = each_pair.rpartition('-')[0] + NORMAL_SIG
            tumor_name = each_pair.rpartition('-')[0] + TUMOR_SIG + '-' + each_pair.rpartition('-')[2]
        else:
            normal_name = each_pair + NORMAL_SIG
            tumor_name = each_pair + TUMOR_SIG
        try:
            germline = glob.glob(os.path.join(csv_output_dir, 'germline', normal_name + '*.csv'))[0]
            somatic_normal = glob.glob(os.path.join(csv_output_dir, 'somatic', normal_name + '*.csv'))[0]
            somatic_tumor = glob.glob(os.path.join(csv_output_dir, 'somatic', tumor_name + '_with_info.csv'))[0]
            somatic_pair = glob.glob(os.path.join(csv_output_dir, 'somatic', each_pair + '_with_info.csv'))[0]
        except:
            import pdb;pdb.set_trace()

    # filter_pipelines2('/home/liaoth/project/170801_XK/result/germline/XK-27W.merged.anno.csv.hg19_multianno.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27W.mt2.merged.anno.hg19_multianno.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27T_with_info.csv',
    #                  '/home/liaoth/project/170801_XK/result/somatic/XK-27_with_info.csv',
    #                  output_file,
    #                  pp=arg_list)