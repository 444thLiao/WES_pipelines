
import sys,os,glob,re
def cal_fastq_reads(fq_path):
    a = os.popen("zgrep '^+$' -c %s" % fq_path)
    return int(a.read())

def cal_fastq_bp(fq_path):
    a = os.popen("zgrep -E '^[ACTGN]+$' $path | wc")
    b = a.read()
    b = b.split(' ')
    return int(b[4])-int(b[1])

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




if __name__ == '__main__':
    setting_file = sys.argv()[-1]
    dir_path = os.path.dirname(setting_file)

    sys.path.insert(0,dir_path)
    from setting import *

    field_names = ['Sample ID','raw data/bp','clean data/bp','genome mapping','target mapping',
                   'remove dup target mapping','dup','target_length/bp','avg_depth','>1X','>5X','>10X','>20X']
    sample_names = [os.path.basename(_i) for _i in glob.glob(os.path.join(base_outpath,'{PN}_result/{PN}*[{n}{t}]'.format(PN=PROJECT_NAME,n=NORMAL_SIG,t=TUMOR_SIG)))]
    raw_path = [[_k for _k in glob.glob(os.path.join(base_inpath, '*%s*R1*.gz') % _i) if filter_str not in _k ][0] for _i in sample_names]
    trim_path = glob.glob(os.path.join(base_outpath, '%s_result/trim_result/*.clean.fq.gz' % PROJECT_NAME))

    bam_path = glob.glob(os.path.join(base_outpath, '%s_result/*/*_sorted.bam' % PROJECT_NAME))

