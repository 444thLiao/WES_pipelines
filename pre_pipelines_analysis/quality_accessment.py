from __future__ import print_function
import sys,os,glob,re
import pandas as pd

Total_pos_num = 65894148


# /home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed

def cal_fastq_reads(fq_path):
    a = os.popen("zgrep '^+$' -c %s" % fq_path)
    return int(a.read())

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

if __name__ == '__main__':
    setting_file = sys.argv[-1]
    dir_path = os.path.dirname(setting_file)
    sys.path.insert(0,dir_path)
    from setting import *

    field_names = ['Sample ID','raw data/bp','clean data/bp','genome mapping','target mapping',
                   'remove dup target mapping','dup','target_length/bp','avg_depth','>1X','>5X','>10X','>20X']
    if not filter_str:
        filter_str = '-999'
    sample_names = [os.path.basename(_i) for _i in glob.glob(
        os.path.join(base_outpath, '{PN}_result/{PN}*[{n}{t}]*'.format(PN=PROJECT_NAME, n=NORMAL_SIG, t=TUMOR_SIG)))]
    raw_path = [[_k for _k in glob.glob(os.path.join(base_inpath, '*%s*R1*.gz' % i ) ) if filter_str not in _k ][0] for i in sample_names]
    trim_path = glob.glob(os.path.join(base_outpath, '%s_result/trim_result/*.clean.fq.gz' % PROJECT_NAME))
    result_df = pd.DataFrame(index=sample_names,columns=field_names[1:])
    # import pdb;pdb.set_trace()
    write_bp(trim_path,'clean data/bp')
    print('completing clean data base pair count.')
    raw_path += [_.replace('R1','R2') for _ in raw_path]
    write_bp(raw_path, 'raw data/bp')
    print('completing raw data base pair count.')
    bam_path = glob.glob(os.path.join(base_outpath, '%s_result/*/*_sorted.bam' % PROJECT_NAME))
    for each in bam_path:
        # tmp = os.popen('/usr/bin/samtools flagstat %s' % each)
        # a = tmp.read()
        # maprate = parse_samtools_info(a,'mapped')
        # result_df.loc[os.path.basename(each).split('_sorted')[0], 'genome mapping'] = float(maprate[0])

        if os.path.isfile(each.replace('_sorted.bam','_sorted_cov.info')):
            #tmp = pd.read_csv(each.replace('_sorted.bam','_sorted_cov.info'),sep='\t',low_memory=False,engine=)
            tmp = open(each.replace('_sorted.bam','_sorted_cov.info')).xreadlines()
            bases = ['A','T','C','G']
            each_pos = []
            count = 0
            for line in tmp:
                if count ==0:
                    idx = [line.split('\t').index(_) for _ in bases]
                else:
                    each_pos.append(sum([int(line.split('\t')[_]) for _ in idx]))
                count += 1
            total_base_count = sum(each_pos)
            result_df.loc[os.path.basename(each).split('_sorted')[0], 'target mapping'] = total_base_count
    print('completing sorted bam(aligned bam) base pair count.')
    bam_path = glob.glob(os.path.join(base_outpath, '%s_result/*/*.recal_reads.bam' % PROJECT_NAME))
    for each in bam_path:
        if os.path.isfile(each.replace('.recal.recal_reads', '_cov.info')):
            tmp = open(each.replace('.recal_reads.bam', '_cov.info')).xreadlines()
            bases = ['A','T','C','G']
            id_labels = ['Chr', 'Posistion']
            each_pos_id = []
            each_pos = []
            count = 0
            for line in tmp:
                if count ==0:
                    idx = [line.split('\t').index(_) for _ in bases]
                    id_id_labels = [line.split('\t').index(_) for _ in id_labels]
                else:
                    pos_id = ';'.join([str(line.split('\t')[_]) for _ in id_id_labels])
                    if pos_id not in each_pos_id:
                        each_pos.append(sum([int(line.split('\t')[_]) for _ in idx]))
                    each_pos_id.append(pos_id)
                count += 1
            total_base_count = sum(each_pos)
            avg_depth = total_base_count / float(len(each_pos_id))
            gt1 = len([_ for _ in each_pos if _ > 1])
            gt5 = len([_ for _ in each_pos if _ > 5])
            gt10 = len([_ for _ in each_pos if _ > 10])
            gt20 = len([_ for _ in each_pos if _ > 20])
            idx_name = os.path.basename(each).split('.dedup')[0]
            result_df.loc[idx_name, 'remove dup target mapping'] = total_base_count
            result_df.loc[idx_name, 'avg_depth'] = avg_depth
            result_df.loc[idx_name, '>1X'] = gt1
            result_df.loc[idx_name, '>5X'] = gt5
            result_df.loc[idx_name, '>10X'] = gt10
            result_df.loc[idx_name, '>20X'] = gt20
    print('completing recal bam(before Calling) base pair count with different depth summary.')
    import pdb;

    pdb.set_trace()
