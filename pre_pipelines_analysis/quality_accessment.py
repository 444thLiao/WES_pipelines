from __future__ import print_function
import sys,os,glob,re
import pandas as pd
import tqdm
import multiprocessing
Total_pos_num = 65894148

# /home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed

def cal_fastq_bp(fq_path):
    a = os.popen("zgrep -E '^[ACTGN]+$' %s | wc" % fq_path)
    b = a.read()
    b = b.split(' ')
    return int(b[2].replace('\n',''))-int(b[1])

def parse_samtools_info(long_str,need ):
    info_list = long_str.split('\n')
    if need == 'mapped':
        parsed = info_list[2]
        result = re.findall(' \(([0-9]{1,2}\.[0-9]{1,2})\%',parsed)
        if result:
            return result
        else:
            raise IOError,parsed
    return None

def write_bp(_path,colname):
    for r1 in _path:
        if 'R1' not in r1:
            continue
        else:
            if r1.count('R1') != 1:
                import pdb;pdb.set_trace()
            r2 = r1.replace('R1','R2')
            r1_bp = cal_fastq_bp(r1)
            r2_bp = cal_fastq_bp(r2)
            sample_n = [_ for _ in sample_names if _ in r1]
            if not sample_n:
                import pdb;pdb.set_trace()
            elif len(sample_n) == 1:
                sample_n = sample_n[0]
            else:
                sample_n = sorted(sample_n)[-1]
            result_df.loc[sample_n,colname] = r1_bp + r2_bp

def sum_bam_bp(each,format,target,depths=False):
    if os.path.isfile(each.partition('.')[0] + '_cov.info'):
        tmp = open(each.partition('.')[0] + '_cov.info').xreadlines()
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
        idx_name = os.path.basename(each).split(format)[0]
        result_df.loc[idx_name, target] = total_base_count
        if depths:
            result_df.loc[idx_name,'>1X'] = len([_ for _ in each_pos if _ > 1])
            result_df.loc[idx_name, '>5X'] = len([_ for _ in each_pos if _ > 5])
            result_df.loc[idx_name, '>10X'] = len([_ for _ in each_pos if _ > 10])
            result_df.loc[idx_name, '>20X'] = len([_ for _ in each_pos if _ > 20])
        if target =='remove dup target mapping':
            result_df.loc[idx_name, 'avg_depth'] = avg_depth

if __name__ == '__main__':
    setting_file = sys.argv[-1]
    dir_path = os.path.dirname(setting_file)
    sys.path.insert(0,dir_path)
    from setting import *

    field_names = ['Sample ID','raw data/bp','clean data/bp','genome mapping','target mapping',
                   'remove dup target mapping','dup','target_length/bp','avg_depth','>1X','>5X','>10X','>20X']
    filter_str = '_somatic'

    sample_names = [os.path.basename(_i) for _i in glob.glob(
        os.path.join(base_outpath, '{PN}_result/{PN}*[{n}{t}]*'.format(PN=PROJECT_NAME, n=NORMAL_SIG, t=TUMOR_SIG))) if filter_str not in _i]
    raw_path = [[_k for _k in glob.glob(os.path.join(base_inpath, '*%s*R1*.gz' % i ) ) if filter_str not in _k ][0] for i in sample_names]
    trim_path = glob.glob(os.path.join(base_outpath, '%s_result/trim_result/*.clean.fq.gz' % PROJECT_NAME))
    result_df = pd.DataFrame(index=sample_names,columns=field_names[1:])
    print('start processing......')
    write_bp(trim_path,'clean data/bp')
    print('completing clean data base pair count.')
    raw_path += [_.replace('R1','R2') for _ in raw_path]
    write_bp(raw_path, 'raw data/bp')
    print('completing raw data base pair count.')
    bam_path = glob.glob(os.path.join(base_outpath, '%s_result/*/*_sorted.bam' % PROJECT_NAME))
    for each in tqdm.tqdm(bam_path):
        sum_bam_bp(each,format='_sorted',target='target mapping')
    print('completing sorted bam(aligned bam) base pair count.')
    bam_path = glob.glob(os.path.join(base_outpath, '%s_result/*/*.recal_reads.bam' % PROJECT_NAME))
    for each in tqdm.tqdm(bam_path):
        print(each)
        sum_bam_bp(each=each, format='.recal',target='remove dup target mapping',depths=True)
    result_df.to_csv(os.path.join(base_outpath,'quality_accessment_raw.csv'))
    print('completing recal bam(before Calling) base pair count with different depth summary.')
    import pdb;pdb.set_trace()
