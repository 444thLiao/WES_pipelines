import sys,os,glob,re
import pandas as pd
from collections import defaultdict
import tqdm
from parse_file_name import fileparser
Total_pos_num = 65894148
# Base on /home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed
# If we change this bed file, Total_pos_num needs to change also.

def cal_fastq_bp(fq_path):
    a = os.popen("zgrep -E '^[ACTGN]+$' %s | wc" % fq_path)
    b = a.read().strip()
    b = [_ for _ in b.split(' ') if _]
    result = int(b[2].replace('\n',''))-int(b[1])
    return result

def parse_samtools_info(long_str,need ):
    info_list = long_str.split('\n')
    if need == 'mapped':
        parsed = info_list[2]
        result = re.findall(' \(([0-9]{1,2}\.[0-9]{1,2})\%',parsed)
        if result:
            return result
        else:
            raise IOError(parsed)
    return None

def sum_bam_bp(each,target,idx_name,depths=False):
    if os.path.isfile(each.partition('.')[0] + '_cov.info'):
        tmp = open(each.partition('.')[0] + '_cov.info').readlines()
        bases = ['A','T','C','G']
        each_pos = []
        count = 0
        for line in tmp:
            if count ==0:
                idx = [line.split('\t').index(_) for _ in bases]
            else:
                each_pos.append(sum([int(line.split('\t')[_]) for _ in idx]))
            count += 1
        print('iterations all row in cov.info')
        total_base_count = sum(each_pos)
        avg_depth = total_base_count / float(count)
        result_df.loc[idx_name, target] = total_base_count
        if depths:
            result_df.loc[idx_name,'>1X'] = len([_ for _ in each_pos if _ > 1])
            result_df.loc[idx_name, '>5X'] = len([_ for _ in each_pos if _ > 5])
            result_df.loc[idx_name, '>10X'] = len([_ for _ in each_pos if _ > 10])
            result_df.loc[idx_name, '>20X'] = len([_ for _ in each_pos if _ > 20])
        if target =='remove dup target mapping':
            result_df.loc[idx_name, 'avg_depth'] = avg_depth

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',dest='input_tab', required=True,type=str)
    parser.add_argument('-o', dest='output_dir', required=True, type=str)
    args = parser.parse_args()

    tab_file = os.path.abspath(args.input_tab)
    odir = os.path.abspath(args.output_dir)
    df = fileparser(tab_file)

    field_names = ['Sample ID',
                   'raw data/bp',
                   'clean data/bp',
                   'target mapping',
                   'remove dup target mapping',
                   'dup',
                   'target_length/bp',
                   'avg_depth',
                   '>1X',
                   '>5X',
                   '>10X',
                   '>20X']
    filter_str = '_somatic'

    sample_names = df.get_attr("sample_name")
    raw_path_R1 = df.get_attr("path_R1")
    raw_path_R2 = df.get_attr("path_R2")
    sample_dict = df.get_output_file_path(odir)

    result_df = pd.DataFrame(index=sample_names, columns=field_names[1:])
    access_files = defaultdict(dict)
    for sid in sample_names:
        access_files[sid]["raw"] = raw_path_R1[sid]
        access_files[sid]["trim"] = sample_dict[sid]["trim_R1"]
        access_files[sid]["sorted_bam"] = sample_dict[sid]["sorted_bam"]
        access_files[sid]["recal_bam"] = sample_dict[sid]["recal_bam"]


    for sid,temp_dict in access_files.items():

        num_raw_reads = cal_fastq_bp(temp_dict["raw"])
        num_trimmed_reads = cal_fastq_bp(temp_dict["trim"])
        sum_bam_bp(temp_dict["sorted_bam"], target='target mapping',
                   idx_name=sid)
        sum_bam_bp(temp_dict["recal_bam"],
                   target='remove dup target mapping',
                   idx_name=sid,
                   depths=True)
        result_df.loc[sid,'raw data/bp'] = num_raw_reads
        result_df.loc[sid, 'clean data/bp'] = num_trimmed_reads

    result_df.loc[:,'target_length/bp'] = Total_pos_num
    result_df.to_csv(os.path.join(odir,'quality_accessment_raw.csv'))
    print('completing recal bam(before Calling) base pair count with different depth summary.')
