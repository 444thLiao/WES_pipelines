import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))

import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from parse_file_name import fileparser

Total_pos_num = 65894148


# check at 20190530
# Base on /home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed
# If we change this bed file, Total_pos_num needs to change also.

def cal_fastq_bp(fq_path):
    a = os.popen("zgrep -E '^[ACTGN]+$' %s | wc" % fq_path)
    b = a.read().strip()
    b = [_ for _ in b.split(' ') if _]
    result = int(b[2].replace('\n', '')) - int(b[1])
    # get total number of character and minus the number of line(\n)
    return result


def parse_samtools_info(long_str, need):
    info_list = long_str.split('\n')
    if need == 'mapped':
        parsed = info_list[2]
        result = re.findall(' \(([0-9]{1,2}\.[0-9]{1,2})\%', parsed)
        if result:
            return result
        else:
            raise IOError(parsed)


def sum_bam_bp(cov_file, target, idx_name, depths=False):
    if not os.path.exists(cov_file):
        raise IOError("cov file isnot generated yet.")
    tmp = open(cov_file).readlines()
    bases = ['A', 'T', 'C', 'G']
    depth_per_pos = []
    query_col_idxs = []
    num_r = 0
    for num_r, line in enumerate(tmp):
        if num_r == 0:
            query_col_idxs = [line.split('\t').index(_)
                              for _ in bases]
        else:
            depth_per_pos.append(sum([int(line.split('\t')[_])
                                      for _ in query_col_idxs]))
    depth_per_pos = np.array(depth_per_pos)
    print('iterations all row in cov.info')
    total_base_count = depth_per_pos.sum()
    avg_depth = total_base_count / float(num_r + 1)
    result_df.loc[idx_name, target] = total_base_count

    if depths:
        result_df.loc[idx_name, '>1X'] = (depth_per_pos > 1).sum()
        result_df.loc[idx_name, '>5X'] = (depth_per_pos > 5).sum()
        result_df.loc[idx_name, '>10X'] = (depth_per_pos > 10).sum()
        result_df.loc[idx_name, '>20X'] = (depth_per_pos > 20).sum()
    if target == 'remove dup target mapping':
        result_df.loc[idx_name, 'avg_depth'] = avg_depth


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_tab', required=True, type=str)
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

    sample_names = df.sid
    raw_path_R1 = df.R1
    raw_path_R2 = df.R2
    sample_dict = df.get_full_info(odir)

    result_df = pd.DataFrame(index=sample_names,
                             columns=field_names[1:])
    access_files = defaultdict(dict)
    for sid in sample_names:
        access_files[sid]["raw"] = raw_path_R1[sid]
        access_files[sid]["trim"] = sample_dict[sid]["trim_R1"]
        access_files[sid]["sorted_cov"] = sample_dict[sid]["sorted_cov"]
        access_files[sid]["recal_cov"] = sample_dict[sid]["recal_cov"]

    for sid, temp_dict in tqdm(access_files.items()):
        # todo: multiprocessing.
        num_raw_reads = cal_fastq_bp(temp_dict["raw"]) *2
        num_trimmed_reads = cal_fastq_bp(temp_dict["trim"]) *2
        sum_bam_bp(temp_dict["sorted_cov"],
                   target='target mapping',
                   idx_name=sid)
        sum_bam_bp(temp_dict["recal_cov"],
                   target='remove dup target mapping',
                   idx_name=sid,
                   depths=True)
        result_df.loc[sid, 'raw data/bp'] = num_raw_reads
        result_df.loc[sid, 'clean data/bp'] = num_trimmed_reads
        result_df.loc[sid, 'dup'] = 1 - result_df.loc[sid,"remove dup target mapping"]/result_df.loc[sid,"target mapping"]
    result_df.loc[:, 'target_length/bp'] = Total_pos_num
    result_df.to_csv(os.path.join(odir,
                                  'quality_accessment_raw.csv'))
    print('completing recal bam(before Calling) base pair count with different depth summary.')
