#################################################################
##############Somatic for fillter script#########################
#################################################################
from pandas import DataFrame as df
import argparse


def single_filter(file_path,
                  output_path='',
                  keep={'Func.refGene': ['exonic', 'exonic;splicing']},
                  drop={'ExonicFunc.refGene': ['synonymous SNV']},
                  other_Info={'6': 'PASS'},
                  disease=True,
                  output_file=True,
                  tlod_filter_complement=False):
    # input after annovar file
    file_df = df.from_csv(file_path, index_col=False)
    #print 'first first:',len(file_df)

    ###keep
    filter_file = df(columns=file_df.columns)
    if keep:
        for i in keep.keys():
            filter_file = file_df[file_df.loc[:, i].isin(keep[i])]
    else:
        filter_file = file_df
    #print 'first second:', len(file_df)
    ####drop
    if drop:
        for i in drop.keys():
            for k in drop[i]:
                filter_file = filter_file[filter_file.loc[:, i] != k]
    #print 'first third:', len(file_df)

    extract_index = []
    process_list = list(filter_file.Otherinfo)
    process_list = [_i.split('\t') for _i in process_list]

    for _idx,i in enumerate(process_list):
        if tlod_filter_complement:
            info_list = i[7].split('=')
            if float(info_list[-1]) >= 4.0:
                extract_index.append(_idx)
        else:
            for key in other_Info.keys():
                if i[int(key)] == other_Info[key]:
                    extract_index.append(_idx)

    filter_file = filter_file.iloc[extract_index, :]
    #print 'first forth:', len(file_df)
    if disease:
        process_list = list(filter_file.Otherinfo)
        process_list = [_i.split('\t') for _i in process_list]
        extract_index = []
        for i in process_list:
            if i[-1].partition(':')[0] == '1/1':
                extract_index.append(process_list.index(i))
        filter_file = filter_file.iloc[extract_index, :]
    if output_file:
        with open(output_path,'w') as f1:
            filter_file.to_csv(f1, index=False)
    else:
        return filter_file


def pair_filter(normal_file,
                tumore_file,
                lowF_file,
                NotInNormal_file,
                keep={'Func.refGene': ['exonic', 'exonic;splicing']},
                drop={'ExonicFunc.refGene': ['synonymous SNV']},
                AQ=30,
                min_F=0.1,
                separte = False):
    # AQ: average quality   major for filter 'not in normal'
    # min_F: minium frequency    major for filter 'low F normal'
    if not separte:
        normal_file_filtered = single_filter(normal_file, disease=False, output_file=False, keep=keep, drop=drop)
        tumore_file_filtered = single_filter(tumore_file, disease=False, output_file=False, keep=keep, drop=drop,
                                             tlod_filter_complement=True)
    else:
        normal_file_filtered = df.from_csv(normal_file, index_col=False)
        tumore_file_filtered = df.from_csv(tumore_file, index_col=False)

    #print 'first:',len(tumore_file_filtered)
    normal_file_filtered_index = ['_'.join([normal_file_filtered.iloc[i, 0], str(int(normal_file_filtered.iloc[i, 1])),
                                            str(int(normal_file_filtered.iloc[i, 2]))]) for i in
                                  range(len(normal_file_filtered))]
    tumore_file_filtered_index = ['_'.join([tumore_file_filtered.iloc[i, 0], str(int(tumore_file_filtered.iloc[i, 1])),
                                            str(int(tumore_file_filtered.iloc[i, 2]))]) for i in
                                  range(len(tumore_file_filtered))]
    #print 'second:', len(tumore_file_filtered_index)
    tumore_file_filtered.index = tumore_file_filtered_index
    normal_file_filtered.index = normal_file_filtered_index

    ########################################
    #not_in_normal = []
    not_in_normal = list(set(tumore_file_filtered_index).difference(set(normal_file_filtered_index)))
    not_in_normal_df = tumore_file_filtered.loc[not_in_normal,:]
    #print 'third: not in normal df', len(not_in_normal_df)

    process_list_T = list(not_in_normal_df.Otherinfo)
    process_list_T = [_i.split('\t') for _i in process_list_T]
    #print 'process_list_T', len(not_in_normal_df)
    '''
    ['chr1',
     '69511',
     'rs75062661',
     'A',
     'G',
     '.',
     'PASS',
     'DB;ECNT=1;HCNT=2;MAX_ED=.;MIN_ED=.;NLOD=0.00;TLOD=51.39',
     'GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1',
     '0/1:0,20:1.00:9:11:0.550:0,625:0:0']
    '''

    T_filtered_low_quality = []
    for i in range(len(process_list_T)):
        cache = dict(zip(process_list_T[i][-2].split(':'), process_list_T[i][-1].split(':')))
        if 'AD' not in cache:
            #print cache, process_list_T[i]
            continue
        if cache['AD']!='.':
            if float(cache['AD'].split(',')[1]) != 0:
                if float(cache['QSS'].split(',')[1]) / float(cache['AD'].split(',')[1]) >= AQ:
                    T_filtered_low_quality.append(i)
        else:
            continue
    #print len(T_filtered_low_quality),len(not_in_normal)
    #print 'T_filtered_low_quality', len(not_in_normal_df)
    not_in_normal_df = not_in_normal_df.iloc[T_filtered_low_quality,]
    #print 'forth: not in normal df', len(not_in_normal_df)
    ###################################################
    both_in_NT = []
    for k in normal_file_filtered_index:
        if k in tumore_file_filtered_index:
            both_in_NT.append(normal_file_filtered_index.index(k))

    normal_file_filtered_first = normal_file_filtered.iloc[both_in_NT,]

    process_list_N = list(normal_file_filtered_first.Otherinfo)
    process_list_N = [_i.split('\t') for _i in process_list_N]

    N_filtered_low_frequency = []
    for i in range(len(process_list_N)):
        cache = dict(zip(process_list_N[i][-2].split(':'), process_list_N[i][-1].split(':')))
        try:
            if float(cache['AD'].split(',')[1]) + float(cache['AD'].split(',')[0]) != 0:
                if float(cache['AD'].split(',')[1]) / (float(cache['AD'].split(',')[1]) + float(cache['AD'].split(',')[0])) <= min_F:
                    N_filtered_low_frequency.append(i)
        except:
            pass

    low_F_normal_df = normal_file_filtered_first.iloc[N_filtered_low_frequency,]

    with open(lowF_file, 'w') as f1:
        low_F_normal_df.to_csv(f1, index=False, sep=',')

    with open(NotInNormal_file, 'w') as f1:
        not_in_normal_df.to_csv(f1, index=False, sep=',')
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-N', '--input1', dest='input_file1', type=str, help='input file1')
    parser.add_argument('-T', '--input2', dest='input_file2', type=str, help='input file2')
    parser.add_argument('-o', '--single_output', dest='single_output_file', type=str, help="output file0")
    parser.add_argument('-lF', '--lowF', dest='lowF_file', type=str, help="output file1")
    parser.add_argument('-nN', '--NotInNormal', dest='NotInNormal_file', type=str, help="output file2")
    parser.add_argument('-oall', '--ALL_in_one', dest='One_file_container', type=str, help="output file3")
    parser.add_argument('-s', '--separte', dest='separte', type=str, help="separte f1 and f1 filter")

    args = parser.parse_args()
    annovar1 = args.input_file1
    annovar2 = args.input_file2

    filter_sep = args.separte

    single_output = args.single_output_file

    lowF_file = args.lowF_file
    NotInNormal_file = args.NotInNormal_file

    outall_on = args.One_file_container
    # print annovar2

    if filter_sep:
        pair_filter(annovar1, annovar2, d, NotInNormal_file, AQ=30,separte = True)
    else:
        if annovar2 == None:
            # print annovar1
            single_filter(annovar1, output_path=single_output, disease=False,keep={})
        elif annovar1 != '' and annovar2 != '':
            pair_filter(annovar1, annovar2, lowF_file, NotInNormal_file,keep={},AQ=30)

    # e.g python somatic_filter_script.py -N '/home/liaoth/project/170123_XK/server_result/somatic/XK-2W_S17.mt2.merged.anno.hg19_multianno.csv' -T '/home/liaoth/project/170123_XK/server_result/somatic/XK-2_merged_TandPair.csv' -lF '/home/liaoth/project/170123_XK/server_result/somatic/filtered/XK-2_mt2_low_F_nofilter.csv' -nN '/home/liaoth/project/170123_XK/server_result/somatic/filtered/XK-2_mt2_NotInNormal_nofilter.csv'