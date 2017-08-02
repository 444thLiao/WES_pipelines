from pandas import DataFrame as df
import os,time
from progressbar import *

"""
Function to add extra column(Majority about mut_percent, ref_cov, mut_cov) from depth_info into csv_result file.



Type: Formatted-for function.
"""



def add_per_info(result_csv,output_csv,tumor_bam,normal_bam):
    '''

    :param result_csv:
    :param output_csv:
    :param tumor_bam:
    :param normal_bam:
    :return:
    '''
    print '{:#^40}'.format('Start Whole project...')
    t1 = time.time()
    normal_cov_info = normal_bam.partition('.')[0]+'_cov.info'
    tumor_cov_info = tumor_bam.partition('.')[0]+'_cov.info'
    ori_csv = df.from_csv(result_csv, index_col=False)
    if not os.path.isfile(normal_cov_info) or not os.path.isfile(tumor_cov_info):
        raise IOError,"cov_info doesn't exist, please use cal Cov pyfile first."
    normal_cov_info_df = df.from_csv(normal_cov_info, index_col=False, sep='\t')
    tumor_cov_info_df = df.from_csv(tumor_cov_info, index_col=False, sep='\t')

    cutted_N_df = normal_cov_info_df[normal_cov_info_df.Posistion.isin(list(ori_csv.Start))]
    cutted_T_df = tumor_cov_info_df[tumor_cov_info_df.Posistion.isin(list(ori_csv.Start))]

    cutted_N_df_index = [';'.join(_pair) for _pair in
                         zip(list(cutted_N_df.Chr), [str(_cache) for _cache in list(cutted_N_df.Posistion)])]
    cutted_T_df_index = [';'.join(_pair) for _pair in
                         zip(list(cutted_T_df.Chr), [str(_cache) for _cache in list(cutted_T_df.Posistion)])]
    cutted_N_df.index = cutted_N_df_index
    cutted_T_df.index = cutted_T_df_index

    added_col = {'N_mut_cov': [],
     'N_mut_per': [],
     'N_ref_cov': [],
     'T_mut_cov': [],
     'T_mut_per': [],
     'T_ref_cov': []}
    t2 = time.time()
    print '{:#^40}'.format('Loaded/Inited all required file...... Using %d ' % (t2-t1))
    print '{:#^40}'.format('Star Iteration.......')
    pbar = ProgressBar().start()
    for _index in list(ori_csv.index):
        Ref = ori_csv.loc[_index, 'Ref']
        Alt = ori_csv.loc[_index, 'Alt']
        Chr = ori_csv.loc[_index, 'Chr']
        Pos = ori_csv.loc[_index, 'Start']
        End = int(ori_csv.loc[_index, 'End'])
        try:
            N_cov_info = cutted_N_df.ix['%s;%i' % (Chr, Pos)]
            T_cov_info = cutted_T_df.ix['%s;%i' % (Chr, Pos)]

        except:
            for _key in added_col:
                added_col[_key].append('Off target')
            continue

        if Ref in ['A', 'C', 'G', 'T'] and Alt in ['A', 'C', 'G', 'T']:
            T_mut_per = T_cov_info.ix[Alt + '_Rate']
            T_ref_cov = T_cov_info.ix[Ref]
            T_mut_cov = T_cov_info.ix[Alt]

            N_mut_per = N_cov_info.ix[Alt + '_Rate']
            N_ref_cov = N_cov_info.ix[Ref]
            N_mut_cov = N_cov_info.ix[Alt]

            added_col['T_mut_per'].append(T_mut_per)
            added_col['T_ref_cov'].append(T_ref_cov)
            added_col['T_mut_cov'].append(T_mut_cov)
            added_col['N_mut_per'].append(N_mut_per)
            added_col['N_ref_cov'].append(N_ref_cov)
            added_col['N_mut_cov'].append(N_mut_cov)

        else:

            exec_pysamstats = '/home/liaoth/tools/pysamstats_venv/bin/pysamstats'
            fasta = '/home/liaoth/data/hg19/ucsc.hg19.fasta'
            t_result = os.popen(
                '{exe} --fasta {fa} --type variation {t_f} -c {chr} -s {start} -e {end} -u'
                .format(exe=exec_pysamstats, fa=fasta, t_f=tumor_bam, chr=Chr, start=Pos, end=End+1)).read()
            n_result = os.popen(
                '{exe} --fasta {fa} --type variation {t_f} -c {chr} -s {start} -e {end} -u'
                .format(exe=exec_pysamstats, fa=fasta, t_f=normal_bam, chr=Chr, start=Pos, end=End+1)).read()

            t_col = t_result.split('\r\n')[0].split('\t')
            t_content = [_row.split('\t') for _row in t_result.split('\r\n')[1:-1]]
            t_data_df = df(columns=t_col, data=t_content,dtype=int)

            n_col = n_result.split('\r\n')[0].split('\t')
            n_content = [_row.split('\t') for _row in n_result.split('\r\n')[1:-1]]
            n_data_df = df(columns=n_col, data=n_content,dtype=int)
            if t_data_df.sum()['reads_pp'] > 0 and len(t_data_df) >0:
                if Ref == '-' or len(Alt) > len(Ref):
                    #insertions
                    added_col['T_mut_per'].append(t_data_df.sum()['insertions_pp']/float(t_data_df.sum()['reads_pp']))
                    added_col['T_ref_cov'].append(t_data_df.sum()['matches_pp']/float(len(t_data_df)))
                    added_col['T_mut_cov'].append(t_data_df.sum()['insertions_pp']/float(len(t_data_df)))
                    added_col['N_mut_per'].append(n_data_df.sum()['insertions_pp']/float(n_data_df.sum()['reads_pp']))
                    added_col['N_ref_cov'].append(n_data_df.sum()['matches_pp']/float(len(n_data_df)))
                    added_col['N_mut_cov'].append(n_data_df.sum()['insertions_pp']/float(len(n_data_df)))

                elif Alt == '-' or len(Ref) > len(Alt):
                    #deletions
                    added_col['T_mut_per'].append(t_data_df.sum()['deletions_pp']/float(t_data_df.sum()['reads_pp']))
                    added_col['T_ref_cov'].append(t_data_df.sum()['matches_pp']/float(len(t_data_df)))
                    added_col['T_mut_cov'].append(t_data_df.sum()['deletions_pp']/float(len(t_data_df)))
                    added_col['N_mut_per'].append(n_data_df.sum()['deletions_pp']/float(n_data_df.sum()['reads_pp']))
                    added_col['N_ref_cov'].append(n_data_df.sum()['matches_pp']/float(len(n_data_df)))
                    added_col['N_mut_cov'].append(n_data_df.sum()['deletions_pp']/float(len(n_data_df)))
                else:
                    for _key in added_col:
                        added_col[_key].append('Wrong pos')
            else:
                for _key in added_col:
                    added_col[_key].append('Off target')
        pbar.update((float(_index)/len(ori_csv))*100)
    pbar.finish()
    for _key in added_col:
        ori_csv.loc[:,_key] = added_col[_key]
    print '{:#^40}'.format('Almost Completing. Iteration used %d.' % (time.time()-t2))
    print '{:#^40}'.format('filtering all unconvinced snp/indel.')
    try:
        a_num = ori_csv[ori_csv.N_mut_per != 'Off target'] # it doesn't inside the bed file.
    except:
        a_num = ori_csv
    a_num.T_mut_per = a_num.T_mut_per.astype(float)
    a_num.N_mut_per = a_num.N_mut_per.astype(float)
    a_num = a_num[a_num.T_mut_per >= a_num.N_mut_per]
    a_num = a_num[a_num.T_mut_per != 0]

    with open(output_csv,'w') as f1:
        a_num.to_csv(f1,index=False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', type=str,required=True,help="input csv file")
    parser.add_argument('-o', '--output', dest='output', type=str, required=True, help="output csv file")
    parser.add_argument('-tb',dest = 'tumor_bam',type=str,required=True,help='path to specify the tumor sample bam file, normally suffix is recal_reads.bam ')
    parser.add_argument('-nb',dest='normal_bam',type=str,required=True,help='path to specify the normal sample bam file, normally suffix is recal_reads.bam')
    args = parser.parse_args()
    input = args.input
    output = args.output
    tb = args.tumor_bam
    nb = args.normal_bam
    add_per_info(input,output,tb,nb)




'''
#normally, a=add_per_info()

#a = a[a.N_mut_per != 'Off target']

#a_indel_index = list(a[a.T_mut_per == 'Indel'].index)
#a_num_index = list(a[a.T_mut_per != 'Indel'].index)
a_num = a
#a_num = a.loc[a_num_index,:]
a_num.T_mut_per = a_num.T_mut_per.astype(float)
a_num.N_mut_per = a_num.N_mut_per.astype(float)

a_num = a_num[a_num.T_mut_per >= a_num.N_mut_per]
a_num = a_num[a_num.T_mut_per != 0]

a_final = a.loc[list(a_num.index),formatt_col]
'''