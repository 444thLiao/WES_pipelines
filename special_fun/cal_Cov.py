import pysam
from pandas import DataFrame as df
import luigi
from progressbar import *
"""
Function to cal each base info from given bed file and bam file.
:type Cal_Function
"""
progress = ProgressBar()
def cal_fun(bam_path, bed_file):

    bed_file = df.from_csv(bed_file, index_col=False, sep='\t')
    try:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
    except:
        print "bamefile need to cal doesn't exist"
        raise IOError

    fastafile = pysam.FastaFile(filename='/home/liaoth/data/hg19/ucsc.hg19.fasta')
    ####define each file path.
    result = 'Gene\tChr\tPosistion\tReference\tbase\tA\tC\tG\tT\tA_Rate\tC_Rate\tG_Rate\tT_Rate\t1\t2\t3\t4\n'
    #define columns.
    print 'STARTING TO ITERATION.'

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
    with open(bam_path.partition('.')[0]+'_cov.info', 'w') as f1:
        f1.write(result)
    print 'Cal cov info complete'


import glob
class cal_Cov(luigi.Task):
    sampleName = luigi.Parameter()
    bed = luigi.Parameter()

    def output(self):
        #sample_name = pfn(self.sampleName, 'sample_name')
        #project_name = pfn(self.sampleName, 'project_name')
        return luigi.LocalTarget(glob.glob('/home/liaoth/project/YZ_XR/output/YZ_XR_result/{_i}*/{_i}*.recal_reads.bam'.format(_i=self.sampleName))[0]
                                 .replace('.recal_reads.bam','_cov.info'))

    def run(self):
        #sample_name = pfn(self.sampleName, 'sample_name')
        #project_name = pfn(self.sampleName, 'project_name')
        #bam_file = output_fmt.format(path=base_outpath, PN=project_name, SN=sample_name) + '.recal_reads.bam'
        bam_file = glob.glob('/home/liaoth/project/YZ_XR/output/YZ_XR_result/{_i}*/{_i}*.recal_reads.bam'.format(_i=self.sampleName))[0]
        cal_fun(bam_file, bed_file=self.bed)


class workflow(luigi.Task):
    x = luigi.Parameter(default='NJZL_S1')
    bed = luigi.Parameter(default='/home/liaoth/project_formal/NYS_miseq/AA_amplicon_info.bed')

    def requires(self):
        samples_IDs = str(self.x).split(',')
        for i in samples_IDs:
            yield cal_Cov(sampleName=i, bed=self.bed)


#from main import *
# python -m luigi --module cal_Cov workflow --x NYS-A_S1,NYS-B_S2,NYS-C_S3,NYS-D_S4 --bed '/home/liaoth/project_formal/170123_NY/AA_amplicon_info.bed' --parallel-scheduling --workers 11 --local-scheduler
