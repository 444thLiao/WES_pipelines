import pysam
from pandas import DataFrame as df
import argparse,time
from progressbar import *
"""
Function to cal each base info from given bed file and bam file.
:type Cal_Function
"""
progress = ProgressBar()
def cal_fun(bam_path, bed_file,REF_file='/home/liaoth/data/hg19/ucsc.hg19.fasta'):
    print 'Start loading required file......'
    t1 = time.time()
    bed_file = df.from_csv(bed_file, index_col=False, sep='\t')
    try:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
    except:
        print "bamefile need to cal doesn't exist"
        raise IOError

    fastafile = pysam.FastaFile(filename=REF_file)
    ####define each file path.
    result = 'Gene\tChr\tPosistion\tReference\tbase\tA\tC\tG\tT\tA_Rate\tC_Rate\tG_Rate\tT_Rate\t1\t2\t3\t4\n'
    #define columns.
    print 'Loading required file. using %d' % (time.time()-t1)
    print 'STARTING TO ITERATION. using '
    t2=time.time()
    total = len(bed_file)
    pbar = ProgressBar().start()
    for i in range(len(bed_file)):
        Chr = bed_file.iloc[i, 0]
        start = min(int(bed_file.iloc[i, 2]), int(bed_file.iloc[i, 1]))
        end = max(int(bed_file.iloc[i, 2]), int(bed_file.iloc[i, 1]))
        Gene_name = str(bed_file.iloc[i, 3])
        #fetch basic info.
        coverage_ACGT = bamfile.count_coverage(Chr, start - 1, end, read_callback='nofilter',
                                               quality_threshold=0)  ###fixed coordinate
        base_counter = dict(zip(['A', 'C', 'G', 'T'], coverage_ACGT))
        ref_base_str = fastafile.fetch(Chr, start - 1, end)
        if sum(sum(j) for j in coverage_ACGT) != 0:
            for base_n in range(start, end + 1):
                n_read = int(sum([k[range(start, end + 1).index(base_n)] for k in coverage_ACGT]))
                ###total base num
                if n_read == 0:
                    continue
                    #if position here didn't have any base, then pass this pos.
                result += Gene_name + '\t' + Chr + '\t' + str(base_n) + '\t' + ref_base_str[
                    range(start, end + 1).index(base_n)] + '\t' + '\t'
                for base in ['A', 'C', 'G', 'T']:
                    result += str(base_counter[base][range(start, end + 1).index(base_n)]) + '\t'
                    #write the A/T/G/C num.

                for base in ['A', 'C', 'G', 'T']:
                    result += str(
                            round(float(base_counter[base][range(start, end + 1).index(base_n)]) / n_read, 4)) + '\t'
                    # wirte the rate of A/T/C/G RATE

                little_rank_list = [(each_one[1][range(start, end + 1).index(base_n)], each_one[0]) for each_one in
                                  base_counter.items()]
                #construct a list to sort A/T/C/G each base num.
                result += '\t'.join([ranked_base[1] for ranked_base in sorted(little_rank_list, reverse=True)]) + '\n'
                #make it columns and tag a '\n' sign.
            if result[-1:] != '\n':
                result += '\n'
            #in case the last base didn't have any base, it will continue, so need to check the last sign.
        else:
            pass
        pbar.update(int((i / (total - 1)) * 100))
    pbar.finish()
    with open(bam_path.partition('.')[0]+'_cov.info', 'w') as f1:
        f1.write(result)
    print 'Cal cov info complete'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', dest='bam_path', type=str,required=True,help="bam file path")
    parser.add_argument('-B','--bed',dest= 'bed_path',type=str,required=True,help='bed file path')
    parser.add_argument('-r', '--ref', dest='ref_fasta', type=str, help='reference fasta file path')
    args = parser.parse_args()

    bam_path = args.bam_path
    bed_path = args.bed_path
    ref_path = args.ref_fasta
    if ref_path:
        cal_fun(bam_path=bam_path,bed_file=bed_path,REF_file=ref_path)
    else:
        cal_fun(bam_path=bam_path, bed_file=bed_path)


