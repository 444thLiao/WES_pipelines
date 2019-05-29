from __future__ import print_function
import pysam
import pandas as pd
import argparse
import tqdm
"""
confirmed 2018.03.19
Function to cal each base info from given bed file and bam file.
:type Cal_Function
"""

def cal_fun(bam_path, bed_file,REF_file='/home/liaoth/data/hg19/ucsc.hg19.fasta'):
    print('Start loading required file......')
    bed_file = pd.read_csv(bed_file, index_col=False, sep='\t', header=None)
    try:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
    except:
        print("bamefile need to cal doesn't exist")
        raise IOError

    fastafile = pysam.FastaFile(filename=REF_file)
    f1 = open(bam_path.partition('.')[0] + '_cov.info', 'w')
    ####define each file path.
    result = 'Gene\tChr\tPosition\tReference\tbase\tA\tC\tG\tT\tA_Rate\tC_Rate\tG_Rate\tT_Rate\t1\t2\t3\t4\n'
    f1.write(result)
    f1.flush()
    #define columns.
    #print 'Loading required file. using %d' % (time.time()-t1)
    print('STARTING TO ITERATION. ')
    count = 0
    pro_count = 0
    if debug_:
        import pdb;pdb.set_trace()
    for i in tqdm.tqdm(range(bed_file.shape[0]),
                       total=bed_file.shape[0]):
        Chr, start, end, Gene_name = bed_file.iloc[i, [0, 1, 2, 3]].values
        start = min(int(start), int(end))
        end = max(int(start), int(end))
        #fetch basic info.
        coverage_ACGT = bamfile.count_coverage(Chr, start, end, read_callback='nofilter',
                                               quality_threshold=0)  ###fixed coordinate
        base_counter = dict(zip(['A', 'C', 'G', 'T'], coverage_ACGT))
        ref_base_str = fastafile.fetch(Chr, start, end)
        if sum(sum(j) for j in coverage_ACGT) != 0:
            for idx_n, base_n in enumerate(range(start, end)):
                pro_count += 1
                n_read = int(sum([k[idx_n] for k in coverage_ACGT]))
                ###total base num
                if n_read == 0:
                    continue
                    #if position here didn't have any base, then pass this pos.
                result = Gene_name + '\t' + Chr + '\t' + str(base_n) + '\t' + ref_base_str[
                    idx_n].upper() + '\t' + '\t'
                for base in ['A', 'C', 'G', 'T']:
                    result += str(base_counter[base][idx_n]) + '\t'
                    count += base_counter[base][idx_n]
                    #write the A/T/G/C num.

                for base in ['A', 'C', 'G', 'T']:
                    result += str(
                        round(float(base_counter[base][idx_n]) / n_read, 4)) + '\t'
                    # wirte the rate of A/T/C/G RATE

                little_rank_list = [(each_one[1][idx_n], each_one[0]) for each_one in
                                    base_counter.items()]
                #construct a list to sort A/T/C/G each base num.
                result += '\t'.join([ranked_base[1] for ranked_base in sorted(little_rank_list, reverse=True)]) + '\n'
                #make it columns and tag a '\n' sign.
                f1.write(result)
                f1.flush()
            #in case the last base didn't have any base, it will continue, so need to check the last sign.
        else:
            pass
    # with open(bam_path.partition('.')[0]+'_cov.info', 'w') as f1:
    #     f1.write(result)
    print('Cal cov info complete.total base is ',count,bam_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', dest='bam_path', type=str,required=True,help="bam file path")
    parser.add_argument('-B','--bed',dest= 'bed_path',type=str,required=True,help='bed file path')
    parser.add_argument('-r', '--ref', dest='ref_fasta', type=str, help='reference fasta file path')
    parser.add_argument('-debug', '--debug', dest='debug_', action="store_true", help="Entering debug mode.")
    args = parser.parse_args()

    bam_path = args.bam_path
    bed_path = args.bed_path
    ref_path = args.ref_fasta
    debug_ = args.debug_
    if ref_path:
        cal_fun(bam_path=bam_path,
                bed_file=bed_path,
                REF_file=ref_path)
    else:
        cal_fun(bam_path=bam_path,
                bed_file=bed_path)


