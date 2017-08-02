import os,pandas,sys
from pandas import DataFrame as df

WES_total = 60700153
WES_BED_PATH = '/home/liaoth/project_formal/170602_XK/WES_sureselect_v6.bed'
def cal_fastq(fastq):
    cache = os.popen("zgrep -c '^+$' %s " % fastq)
    bp_count = sys.system("zcat {input} | paste - - - - | cut -f 2 | tr -d '\n' | wc - c")
    return int(cache.read())


def cal_sam(bam):
    if os.path.isfile(bam.partition('.')[0]+'_cov.info'):
        pass
    else:
        os.system('python /home/liaoth/tools/Whole_pipelines/special_fun/cal_Cov_script_version.py -b %s -B %s'.format(
            bam,WES_BED_PATH))
    raw_info = pandas.read_csv(bam.partition('.')[0] + '_cov.info', sep='\t', index_col=False, engine='python')
    raw_info.base = list(raw_info.loc[:, ['A', 'T', 'C', 'G']].sum(1))
    current_bp = raw_info.sum(0).ix['base']
    return current_bp

def series_output(default=['raw_data','clean_data','mapping_rate','target_mapping','dedup_target_mapping',]):
    pass