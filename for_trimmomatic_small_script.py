import os
import glob
#example:
#input1 = 'XK-9F_S9_L001_R1_001'
#input2 = 'XK-9F_S9_L001_R2_001'
#base_inpath = '/media/vd1/170519_fqs
#base_outpath = '/home/liaoth/project_formal/170159_fqs/trimed_result/Q30'

if __name__ == '__main__':
    path_data_dir = '/home/liaoth/data/ck/raw'
# = raw_input("Entering the path of data file you store. e.g:/home/liaoth/project/mbiom/raw_datas.   :")
    path_output_dir = '/home/liaoth/data/ck/trimmed_result/'
# = raw_input("Entering the path of the output data you want to store. e.g:/home/liaoth/project/mbiom/trimmed_result.   ")
    if not os.path.isdir(path_data_dir):
        os.makedirs(path_data_dir)
    if not os.path.isdir(path_output_dir):
        os.makedirs(path_output_dir)

    sample_ids = set([path.rpartition('/')[2].rpartition('_L001_R')[0] for path in glob.glob(path_data_dir+'/*.gz')])
    base_inpath = path_data_dir
    for sample_id in sample_ids:
        input1 = sample_id+'_L001_R1_001'
        input2 = sample_id+'_L001_R2_001'

        os.system("java -jar /home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 5 {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:100".format(
                    input1=input1, input2=input2, base_in=base_inpath, base_out=path_output_dir,
                    output=path_output_dir + input1.replace('_L001_R1_001','_trimed.log')))

for each in `ls /home/liaoth/project/XK_WES/180309_all/ output / XK_result / * / *.recal_reads.bam`; do python ~ / tools / Whole_pipelines / pre_pipelines_analysis / cal_Cov_script_version.py -b $each -B / home / liaoth / data_bank / XK_WES / Sureselect_V6_COSMIC_formal.bed -r / home / db_public / hg19 / ucsc.hg19.fasta & done
