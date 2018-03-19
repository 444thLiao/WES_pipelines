import pandas as pd
import tqdm
import time,os
import pysamstats
import pysam
"""
Function to add extra column(Majority about mut_percent, ref_cov, mut_cov) from depth_info into csv_result file.



Type: Formatted-for function.
"""

tmp = ['deletions',
 'reads_pp',
 'reads_all',
 'matches_pp',
 'pos',
 'insertions_pp',
 'mismatches',
 'C_pp',
 'matches',
 'N_pp',
 'ref',
 'T_pp',
 'A',
 'C',
 'G',
 'deletions_pp',
 'N',
 'mismatches_pp',
 'T',
 'A_pp',
 'chrom',
 'insertions',
 'G_pp']

def parse_bam(bam,Chr,Pos,End,Ref,Alt,REF,sig='N',):
    #exec_pysamstats = '/home/liaoth/tools/pysamstats_venv/bin/pysamstats'
    #
    # t_result = os.popen(
    #     '{exe} --fasta {fa} --type variation {t_f} -c {chr} -s {start} -e {end} -u'
    #         .format(exe=exec_pysamstats, fa=fasta, t_f=bam_path, chr=Chr, start=Pos, end=End + 1)).read()
    result = pysamstats.stat_variation(bam,REF,chrom=Chr,start=Pos,end=End+1,truncate=True,one_based=True)
    result_dict = {}
    for idx,record in enumerate(result):
        result_dict[idx] = record
    t_data_df = pd.DataFrame.from_dict(result_dict,orient='index')

    if t_data_df.shape[0] == 0:
        t_data_df = t_data_df.append(pd.DataFrame(columns=tmp))
        t_data_df.loc[0, :] = 0

    if Ref in ['A', 'C', 'G', 'T'] and Alt in ['A', 'C', 'G', 'T']:
        T_ref_cov,T_mut_cov,T_cov = t_data_df.loc[0, [Ref + '_pp',Alt + '_pp','reads_pp']]

        if T_cov != 0:
            T_mut_per = T_mut_cov / float(T_cov)
        else:
            T_mut_per = 0

        added_col['%s_mut_per' % sig].append(T_mut_per)
        added_col['%s_ref_cov' % sig].append(T_ref_cov)
        added_col['%s_mut_cov' % sig].append(T_mut_cov)
    else:
        if t_data_df.sum()['reads_pp'] > 0 and len(t_data_df) > 0:
            if Ref == '-' or len(Alt) > len(Ref):
                # insertions
                added_col['%s_mut_per' % sig].append(t_data_df.sum()['insertions_pp'] / float(t_data_df.sum()['reads_pp']))
                added_col['%s_ref_cov' % sig].append(t_data_df.sum()['matches_pp'] / float(len(t_data_df)))
                added_col['%s_mut_cov' % sig].append(t_data_df.sum()['insertions_pp'] / float(len(t_data_df)))
            elif Alt == '-' or len(Ref) > len(Alt):
                # deletions
                added_col['%s_mut_per' % sig].append(t_data_df.sum(0)['deletions_pp'] / float(t_data_df.sum()['reads_pp']))
                added_col['%s_ref_cov' % sig].append(t_data_df.sum(0)['matches_pp'] / float(len(t_data_df)))
                added_col['%s_mut_cov' % sig].append(t_data_df.sum(0)['deletions_pp'] / float(len(t_data_df)))
            else:
                for _key in [_ for _ in added_col.keys() if sig in _]:
                    added_col[_key].append('Wrong pos')
        else:
            for _key in [_ for _ in added_col.keys() if sig in _]:
                added_col[_key].append('Off target')


def add_per_info(result_csv,output_csv,tumor_bam,normal_bam):
    '''
    :param result_csv:
    :param output_csv:
    :param tumor_bam:
    :param normal_bam:
    :return:
    '''
    fasta = '/home/db_public/hg19/ucsc.hg19.fasta'
    print '{:#^40}'.format('Start Whole project...')
    t1 = time.time()
    tb_info =
    nb_info =
    tb = pysam.AlignmentFile(tumor_bam)
    nb = pysam.AlignmentFile(normal_bam)
    ref_ = pysam.FastaFile(fasta)

    ori_csv = pd.read_csv(result_csv, index_col=None)
    t2 = time.time()
    print '{:#^40}'.format('Loaded/Inited all required file...... Using %d ' % (t2 - t1))
    print '{:#^40}'.format('Star Iteration.......')
    added_col = {'N_mut_cov': [],
     'N_mut_per': [],
     'N_ref_cov': [],
     'T_mut_cov': [],
     'T_mut_per': [],
     'T_ref_cov': []}
    global added_col
    for _index in tqdm.tqdm(list(ori_csv.index)):
        Ref = ori_csv.loc[_index, 'Ref']
        Alt = ori_csv.loc[_index, 'Alt']
        Chr = ori_csv.loc[_index, 'Chr']
        Pos = ori_csv.loc[_index, 'Start']
        End = int(ori_csv.loc[_index, 'End'])
        if len(set(range(Pos,End+1)).difference(range_list)) == End - Pos +1:
            for _key in added_col:
                added_col[_key].append('Totally Off target')
            continue
        if tumor_bam:
            parse_bam(tb,Chr,Pos,End,Ref,Alt,ref_,sig='T')
        if normal_bam:
            parse_bam(nb, Chr, Pos, End, Ref, Alt,ref_, sig='N')
    for _key in added_col:
        ori_csv.loc[:,_key] = added_col[_key]
    print '{:#^40}'.format('Almost Completing. Iteration used %d.' % (time.time()-t2))
    print '{:#^40}'.format('filtering all unconvinced snp/indel.')
    try:
        a_num = ori_csv[ori_csv.N_mut_per != 'Off target'] # it doesn't inside the bed file.
    except:
        a_num = ori_csv
    # a_num.T_mut_per = a_num.T_mut_per.astype(float)
    # a_num.N_mut_per = a_num.N_mut_per.astype(float)
    # a_num = a_num[a_num.T_mut_per >= a_num.N_mut_per]
    # a_num = a_num[a_num.T_mut_per != 0]

    with open(output_csv,'w') as f1:
        a_num.to_csv(f1,index=False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', type=str,required=True,help="input csv file")
    parser.add_argument('-o', '--output', dest='output', type=str, required=True, help="output csv file")
    parser.add_argument('-tb',dest = 'tumor_bam',type=str,help='path to specify the tumor sample bam file, normally suffix is recal_reads.bam')
    parser.add_argument('-nb',dest='normal_bam',type=str,help='path to specify the normal sample bam file, normally suffix is recal_reads.bam')
    parser.add_argument('-b',dest='bed_file',type=str,required=True,help='path to specify the bed file')
    args = parser.parse_args()
    input = args.input
    output = args.output
    tb = args.tumor_bam
    nb = args.normal_bam
    add_per_info(input,output,tb,nb)