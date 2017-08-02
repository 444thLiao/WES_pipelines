from collections import defaultdict
import luigi
import json

def multi_threads_mt2_single(INBAM, N):
    """
    Receive a bam file to analyze the reads number in bam, then use it to generate N intervals files.
    :param INBAM:
    :param N:
    :return: a list with N tuple, each tuple [0] is id number, tuple[1] is another tuple with multi chrom.
    """
    ###############   make chr list, extract different combined chr in order to construct a small bam.###########################
    bam_stats = os.popen('''/usr/bin/samtools idxstats {bam_file} | cut -f 3'''.format(bam_file=INBAM))
    bam_chr = os.popen('''/usr/bin/samtools idxstats {bam_file} | cut -f 1'''.format(bam_file=INBAM))
    bam_stats = bam_stats.read().splitlines()
    bam_stats = [int(i) for i in bam_stats][:-1]
    bam_chr = bam_chr.read().splitlines()

    each_chr_reads = sorted([list(i) for i in zip(bam_stats, bam_chr)], reverse=True)
    total_reads = sum([int(i) for i in bam_stats])

    each_bucket_size = total_reads / N

    bam_spliting_dict = defaultdict(list)

    index_for_dict = 1
    while len(bam_spliting_dict) < N and len(each_chr_reads) != 0:
        if len(each_chr_reads) == 0:
            break
        biggest_one = each_chr_reads[0]
        each_chr_reads = each_chr_reads[1:]
        if biggest_one[0] >= each_bucket_size:
            bam_spliting_dict[str(index_for_dict)] = [biggest_one[1]]
        else:
            cache_num = biggest_one[0]
            cache_bucket = [biggest_one[1]]
            while cache_num <= each_bucket_size and len(each_chr_reads) != 0:
                smallest_one = each_chr_reads[-1]
                cache_num += smallest_one[0]
                cache_bucket += [smallest_one[1]]
                each_chr_reads = each_chr_reads[:-1]
            bam_spliting_dict[str(index_for_dict)] = cache_bucket
        index_for_dict += 1
    for i in bam_spliting_dict:
        bam_spliting_dict[i] = tuple(bam_spliting_dict[i])
    return bam_spliting_dict.items()


###########BaseRecalibrator and HaplotypeCaller

class QC_trimmomatic(luigi.Task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)

    def output(self):
        sample_name = pfn(self.PE1, 'sample_name')
        project_name = pfn(self.PE1, 'project_name')
        return luigi.LocalTarget(
            '{base}/{PN}_result/trim_result/{SN}_trimed.log'.format(base=base_outpath, PN=project_name, SN=sample_name))

    def run(self):
        project_name = pfn(self.PE1, 'project_name')

        if os.path.isdir('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name)) != True:
            os.makedirs('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name))

        input1 = self.PE1
        input2 = self.PE2
        if input2 != None:
            os.system(
                "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                    input1=input1, input2=input2, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                    output=self.output().path))
        else:
            os.system(
                "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 {base_in}/{input1}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                    input1=input1, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                    output=self.output().path))


class GenerateSam_pair(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        if Pair_data:
            input1 = PE1_fmt.format(input=self.sampleID)
            input2 = PE2_fmt.format(input=self.sampleID)
            return QC_trimmomatic(PE1=input1, PE2=input2)
        else:
            input1 = self.sampleID
            return QC_trimmomatic(PE1=input1)

    def output(self):
        sample_name = self.sampleID
        project_name = pfn(self.sampleID, 'project_name')

        return luigi.LocalTarget(
            output_fmt.format(path=base_outpath, PN=project_name, SN=sample_name) + '.sam')

    def run(self):
        sample_name = self.sampleID
        project_name = pfn(self.sampleID, 'project_name')

        if 'pair_ID' in pfn(self.sampleID, 'all'):
            input1 = '{base_in}/{pe1_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe1_fmt=PE1_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            input2 = '{base_in}/{pe2_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe2_fmt=PE2_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            if not os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) == True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            os.system(
                "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {o}".format(
                    SN=sample_name, REF=REF_file_path, i1=input1, i2=input2, o=self.output().path))
        else:
            input1 = '{base_in}/{SE_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                             SE_fmt=SE_fmt.format(input=self.sampleID),
                                                             input=self.sampleID)
            if os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) != True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            os.system(
                "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1}  > {o}".format(
                    SN=sample_name, REF=REF_file_path, i1=input1, o=self.output().path))


class Convertbam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [GenerateSam_pair(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path[:-3] + 'bam')

    def run(self):
        os.system("samtools view -F 0x100 -bSu %s -o %s" % (self.input()[0].path, self.output().path))


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Convertbam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path[:-4] + '_sorted.bam')

    def run(self):
        os.system("samtools sort -m 100G -f -@ 30 %s %s" % (self.input()[0].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


#########2
class MarkDuplicate(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path[:-11] + '.dedup.bam')

    def run(self):
        if PCR_ON:
            pass
        else:
            os.system(
                "java -Xmx2g -jar ~/tools/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true" % (
                    self.input()[0].path, self.output().path, self.output().path.rpartition('/')[0]))


#########3
class RealignerTargetCreator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.dedup.bam')[0] + '.realign.intervals')

    def run(self):
        os.system(
            "java -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R %s -I %s --known %s -o %s" % (
                REF_file_path, self.input()[0].path, known_gold_cvf, self.output().path))


#########4
class IndelRealigner(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID), RealignerTargetCreator(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.dedup.bam')[0] + '.realign.bam')

    def run(self):
        os.system(
            "java -Xmx5g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.realign.bam')[0] + '.recal_data.table')

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 25 -R %s -I %s -knownSites %s -knownSites %s -o %s" % (
                REF_file_path, self.input()[0].path, db_snp, known_gold_cvf, self.output().path))


#########6
class PrintReads(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID), BaseRecalibrator(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.realign.bam')[0] + '.recal_reads.bam')

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads -R %s -I %s -BQSR %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


########################tumor and normal merged step


#########somatic pipeline

class MuTect2_pair(luigi.Task):
    sample_IDs = luigi.Parameter()

    ####LCB
    def requires(self):
        sampleIDs = self.sample_IDs.split(',')
        return [PrintReads(sampleID=sampleIDs[0]), PrintReads(sampleID=sampleIDs[1])]

    def output(self):
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath, PN=Project_ID, PairN=pair_name) + '.mt2_finished'
        return luigi.LocalTarget(output_path)

    def run(self):
        multi_list = multi_threads_mt2_single(self.input()[0].path, 10)
        #######random choose , here choose N's intervals


        yield [MuTect2_pair_multi_threads(bam_path_N=self.input()[0].path, bam_path_T=self.input()[1].path,
                                          all_info=json.dumps(i), sample_IDs=self.sample_IDs) for i in multi_list]

        os.system('touch {file}'.format(file = self.output().path))


class MuTect2_single(luigi.Task):
    sample_NT = luigi.Parameter()

    ####LCN/TB
    def requires(self):
        return [PrintReads(sampleID=self.sample_NT)]

    def output(self):

        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')
        output_path = somatic_single_output_fmt.format(path=base_outpath, PN=Project_ID, SN=sample_name) + '.mt2_finished'
        return luigi.LocalTarget(output_path)

    def run(self):
        multi_list = multi_threads_mt2_single(self.input()[0].path, 10)
        yield [MuTect2_single_multi_threads(bam_path=self.input()[0].path, all_info=json.dumps(i),
                                            sample_NT=self.sample_NT) for i in multi_list]

        os.system('touch {file}'.format(file = self.output().path))


class Merge_multi_vcf(luigi.Task):
    sample_NT = luigi.Parameter()

    def requires(self):
        if self.sample_NT in pair_bucket:
            pair_value = pair_bucket[self.sample_NT]
            return [MuTect2_pair(sample_IDs=','.join(pair_value)),PrintReads(sampleID=pair_value[0])]
        # doesn't know the tumor/normal order , but it doesn't care.
        else:
            return [MuTect2_single(sample_NT=self.sample_NT),PrintReads(sampleID=self.sample_NT)]

    def output(self):
        output_path = self.input()[0].path.replace('.mt2_finished','_mt2.vcf')
        return luigi.LocalTarget(output_path)

    def run(self):
        input_path = self.input()[0].path

        base_vcf_path = input_path[:input_path.rindex('.')]
        multi_list = multi_threads_mt2_single(self.input()[1].path, 10)

        cache_merge_vcf = [base_vcf_path +'.'+ str(i[0]) for i in multi_list]
        if not self.sample_NT in pair_bucket:
            cache_merge_vcf = [i + '_mt2.vcf' for i in cache_merge_vcf]
        else:
            cache_merge_vcf = [i + '_mt2_combined.vcf' for i in cache_merge_vcf]

        for each_vcf in cache_merge_vcf:
            os.system('/home/liaoth/tools/bcftools/bcftools view {i_vcf} -o {i_vcf}.gz -O z'.format(i_vcf=each_vcf))
            os.system('/home/liaoth/tools/bcftools/bcftools index {i_vcf} -o {i_vcf}.gz -O z'.format(i_vcf=each_vcf))
        os.system('/home/liaoth/tools/bcftools/bcftools concat {vcf_list} -o {final_vcf}'.format(vcf_list='.gz '.join(cache_merge_vcf),
                                                                                    final_vcf=self.output().path))

class Annovar1(luigi.Task):
    sample_ID = luigi.Parameter()

    ####LC(N/T)B
    def requires(self):
        return [Merge_multi_vcf(sample_NT=self.sample_ID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.vcf')[0] + '.merged.av')

    def run(self):
        # prefix=self.input().path.rpartition('_mt2_combined.bam')[0]
        os.system("convert2annovar.pl %s --includeinfo -format vcf4 > %s" % (
            self.input()[0].path, self.output().path))


class Annovar2(luigi.Task):
    sample_ID = luigi.Parameter()

    ####LCB
    def requires(self):
        return [Annovar1(sample_ID=self.sample_ID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.merged.av')[0] + '.merged.anno')

    def run(self):
        os.system(
            "table_annovar.pl %s ~/tools/annovar/humandb/ -buildver hg19 --remove --otherinfo -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,exac03,snp138,ljb26_all,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f -nastring . --otherinfo --csvout --thread 10 --outfile %s --argument '-exonicsplicing -splicing 25',,,,,,,," % (
                self.input()[0].path, self.input()[0].path.rpartition('.merged.av')[0] + '.merged.anno'))
        # os.system("table_annovar.pl %s ~/tools/annovar/humandb/ -buildver hg19 --remove --otherinfo -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,exac03,snp138,ljb26_all,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f -nastring . --otherinfo --csvout --thread 10 --outfile %s --argument '-exonicsplicing -splicing 25',,,,,,,," % (self.input()[2].path,self.input()[2].path.rpartition('.merged.av')[0]+'.merged.anno'))


class workflow(luigi.Task):
    x = luigi.Parameter()

    def requires(self):
        samples_IDs = str(self.x).split(',')

        pair_bucket = defaultdict(list)
        for _x in samples_IDs:
            pair_bucket[pfn(_x, 'pair_name')].append(_x)
        global pair_bucket
        ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}

        samples_IDs += [_x for _x in pair_bucket.keys()]
        for i in samples_IDs:
            yield Annovar2(sample_ID=i)

if __name__ == '__main__':
    luigi.run()


# python -m luigi --module Somatic_multi_threads_mt2 workflow --x LCB --parallel-scheduling --workers 8 --local-scheduler








###############multi mt2 single part##################

class generate_bam(luigi.Task):
    all_info = luigi.Parameter()
    bam_path = luigi.Parameter()
    sample_NT = luigi.Parameter()

    def output(self):
        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')

        output_path = somatic_single_output_fmt.format(path=base_outpath, PN=Project_ID, SN=sample_name)
        intervals_path = output_path + '.' + \
                         json.loads(self.all_info)[0] + '.intervals'
        ####new_dir+'.'+1/2/3/4/5...+'.intervals'
        return luigi.LocalTarget(intervals_path)

    def run(self):
        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')

        output_path = somatic_single_output_fmt.format(path=base_outpath, PN=Project_ID, SN=sample_name)
        intervals_path = output_path + '.' + \
                         json.loads(self.all_info)[0] + '.intervals'

        if not os.path.isdir(output_path.rpartition('/')[0]):
            os.makedirs(output_path.rpartition('/')[0])
        with open(intervals_path, 'w') as f1:
            f1.write('\n'.join(json.loads(self.all_info)[1]))

class MuTect2_single_multi_threads(luigi.Task):
    all_info = luigi.Parameter()
    bam_path = luigi.Parameter()
    sample_NT = luigi.Parameter()

    def requires(self):
        return generate_bam(bam_path=self.bam_path, all_info=self.all_info, sample_NT=self.sample_NT)

    def output(self):
        prefix = self.input().path.replace('.intervals','')

        return luigi.LocalTarget('{prefix}_mt2.vcf'.format(prefix=prefix))

    def run(self):
        prefix = self.input().path.rpartition('.intervals')[0]
        ####~~~~/LCNB.recal_reads.5.bam
        os.system(
            '''java -Xmx10g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T MuTect2 -L {region} --forceActive -out_mode EMIT_ALL_SITES --allSitePLs --artifact_detection_mode --dontTrimActiveRegions --disableOptimizations -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_tumor} --out {prefix}_mt2.vcf --bamOutput {prefix}_mt2.bam --log_to_file {prefix}_mt2.log '''.format(
                region=self.input().path, REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_tumor=self.bam_path,
                prefix=prefix))


class MuTect2_pair_multi_threads(luigi.Task):
    all_info = luigi.Parameter()
    bam_path_N = luigi.Parameter()
    bam_path_T = luigi.Parameter()
    sample_IDs = luigi.Parameter()

    def requires(self):
        sampleIDs = self.sample_IDs.split(',')
        if pfn(sampleIDs[0], 'mt2_for') == NORMAL_SIG:
            normal = sampleIDs[0]
            tumor = sampleIDs[1]
        elif pfn(sampleIDs[0], 'mt2_for') == TUMOR_SIG:
            normal = sampleIDs[1]
            tumor = sampleIDs[0]
        else:
            tumor = ''
            normal = ''
        return [generate_bam(bam_path=self.bam_path_N, all_info=self.all_info, sample_NT=normal),
                generate_bam(bam_path=self.bam_path_T, all_info=self.all_info, sample_NT=tumor)]

    def output(self):
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath, PN=Project_ID, PairN=pair_name)

        ###extract path without N/T
        return luigi.LocalTarget('{prefix}.{ID}_mt2_combined.vcf'.format(prefix=output_path,ID=json.loads(self.all_info)[0]))

    def run(self):
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath, PN=Project_ID, PairN=pair_name)
        intervals_path = self.input()[0].path
        ####random choose , here choose N's intervals

        os.system(
            "java -Xmx10g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T MuTect2 -L {region} --forceActive -out_mode EMIT_ALL_SITES --allSitePLs --dontTrimActiveRegions --disableOptimizations -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:normal {input_normal} --input_file:tumor {input_tumor} --out {prefix}_{ID}_mt2_combined.vcf --bamOutput {prefix}_{ID}_mt2_combined.bam --log_to_file {prefix}_{ID}_mt2_combined.log ".format(
                region=intervals_path, REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_normal=self.bam_path_N,
                input_tumor=self.bam_path_T, prefix=output_path,ID=json.loads(self.all_info)[0]))
        # os.system('touch %s ' % self.output().path)


# class multi_mt2_workflow(luigi.Task):
#     input_path = luigi.Parameter()
#
#     def requires(self):
#         multi_list = multi_threads_mt2_single(self.input_path, 10)
#         for i in multi_list:
#             yield MuTect2_single_multi_threads(bam_path=self.input_path, all_info=i)
#
#     def output(self):
#         return luigi.LocalTarget(
#             '{file}.multi_threads_mt2_finished_1'.format(file=self.input()[0].path.rpartition('/')[0]))
#
#     def run(self):
#         os.system('touch {file}.multi_threads_mt2_finished_1'.format(file=self.input_path.rpartition('/')[0]))

#
from main import *