###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample 
### @GPZ-bioinfo, 20170301
###############################################################################################

import luigi
import glob, time
from main import *


def record_cmdline(message, default=base_outpath + '/germline_pipelines.log'):
    if os.path.isfile(default):
        with open(default, 'a') as f1:
            f1.write(time.ctime() + ' ' * 4 + message + '\n')
    else:
        with open(default, 'w') as f1:
            f1.write('{:#^40}'.format('Starting the somatic pipelines.'))
            f1.write(time.ctime() + ' ' * 4 + message + '\n')


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
            cmdline = "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                input1=input1, input2=input2, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                output=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            cmdline = "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 {base_in}/{input1}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                input1=input1, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                output=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)


class GenerateSam_pair(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        if Pair_data:
            if not self_adjust_fn:
                input1 = PE1_fmt.format(input=self.sampleID)
                input2 = PE2_fmt.format(input=self.sampleID)
            else:
                input_list = glob.glob(base_inpath + '/*' + self.sampleID + '*')
                if filter_str:
                    input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
                input1 = [_i for _i in input_list if R1_INDICATOR in _i][0]
                input2 = [_i for _i in input_list if R2_INDICATOR in _i][0]
            return QC_trimmomatic(PE1=os.path.basename(input1), PE2=os.path.basename(input2))
        else:
            if not self_adjust_fn:
                input1 = SE_fmt.format(input=self.sampleID)
            else:
                input_list = glob.glob(base_inpath + '/*' + self.sampleID + '*')
                if filter_str:
                    input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
                input1 = [_i for _i in input_list if R1_INDICATOR in _i][0]
            return QC_trimmomatic(PE1=os.path.basename(input1))

    def output(self):
        sample_name = pfn(self.sampleID, 'sample_name')
        project_name = pfn(self.sampleID, 'project_name')

        return luigi.LocalTarget(
            output_fmt.format(path=base_outpath, PN=project_name, SN=sample_name) + '.sam')

    def run(self):
        sample_name = pfn(self.sampleID, 'sample_name')
        project_name = pfn(self.sampleID, 'project_name')

        if Pair_data:
            input1 = '{base_in}/{pe1_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe1_fmt=PE1_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            input2 = '{base_in}/{pe2_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe2_fmt=PE2_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            if not os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) == True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {o}".format(
                SN=sample_name, REF=REF_file_path, i1=input1, i2=input2, o=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            input1 = '{base_in}/{SE_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                             SE_fmt=SE_fmt.format(input=self.sampleID),
                                                             input=self.sampleID)
            if os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) != True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1}  > {o}".format(
                SN=sample_name, REF=REF_file_path, i1=input1, o=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)


class Convertbam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [GenerateSam_pair(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.sam', '.bam'))

    def run(self):
        cmdline = "samtools view -F 0x100 -bSu %s -o %s" % (self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Convertbam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.bam', '_sorted.bam'))

    def run(self):
        cmdline = "samtools sort -m 60G -f -@ 30 %s %s" % (self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)
        cmdline = 'samtools index %s' % self.output().path
        os.system(cmdline)
        record_cmdline(cmdline)


#########2
class MarkDuplicate(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        if PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "gatk MarkDuplicates --java-options '-Xmx4g' --INPUT %s --OUTPUT %s --METRICS_FILE %s/dedup_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true -AS true" % (self.input()[0].path, self.output().path, self.output().path.rpartition('/')[0])
        os.system(cmdline)
        record_cmdline(cmdline)

#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam', '.recal_data.table'))

    def run(self):
        cmdline = "gatk BaseRecalibrator --java-options '-Xmx4g' --reference %s --input %s --known-sites %s --known-sites %s --output %s" % (
            REF_file_path, self.input()[0].path, db_snp, known_gold_cvf, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########6
class PrintReads(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID), BaseRecalibrator(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.recal_reads.bam'))

    def run(self):
        cmdline = "gatk ApplyBQSR --java-options '-Xmx4g' --reference %s --input %s --bqsr-recal-file %s --output %s" % (
            REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)
        cmdline = 'samtools index %s' % self.output().path
        os.system(cmdline)
        record_cmdline(cmdline)


#########7
class HaplotypeCaller(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [PrintReads(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.recal_reads.bam', '.raw_variants.vcf'))

    def run(self):
        if bed_file_path != '':

            cmdline = "gatk HaplotypeCaller --java-options '-Xmx4g' --native-pair-hmm-threads 30 --reference {ref} --input {input} --genotyping-mode DISCOVERY --dbsnp {dbsnp} -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {output} --intervals {tar_bed}".format(ref=REF_file_path,input=self.input()[0].path,dbsnp=db_snp,output=self.output().path,tar_bed=bed_file_path)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            cmdline = "gatk HaplotypeCaller --java-options '-Xmx4g' --native-pair-hmm-threads 30 --reference {ref} --input {input} --genotyping-mode DISCOVERY --dbsnp {dbsnp} -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {output}".format(ref=REF_file_path,input=self.input()[0].path,dbsnp=db_snp,output=self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########9
class SelectVariants_a(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [HaplotypeCaller(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_variants.vcf', '.raw_snps.vcf'))

    def run(self):
        cmdline = "gatk SelectVariants --java-options '-Xmx4g' -R %s -V %s -select-type SNP -O %s" % (
            REF_file_path, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########10
class VariantFiltration_a(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [SelectVariants_a(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_snps.vcf', '.filter_snps.vcf'))

    def run(self):
        cmdline = "gatk VariantFiltration --java-options '-Xmx4g' -R %s -V %s --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'my_snp_filter' -O %s" % (
            REF_file_path, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########11
class SelectVariants_b(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [HaplotypeCaller(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_variants.vcf', '.raw_indels.vcf'))

    def run(self):
        cmdline = "gatk SelectVariants --java-options '-Xmx4g' -R %s -V %s -select-type INDEL -O %s" % (
            REF_file_path, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########12
class VariantFiltration_b(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [SelectVariants_b(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_indels.vcf', '.filter_indels.vcf'))

    def run(self):
        cmdline = "gatk VariantFiltration --java-options '-Xmx4g' -R {ref} -V {input} --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filter-name 'my_indel_filter' -O {output}".format(ref=REF_file_path,input=self.input()[0].path,output=self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########13
class CombineVariants(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [VariantFiltration_a(sampleID=self.sampleID), VariantFiltration_b(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace('.filter_indels.vcf', '.merged.vcf'))

    def run(self):
        cmdline = """gatk MergeVcfs --java-options "-Xmx4g"  -R {ref} --INPUT {input_indel} --INPUT {input_snps} --OUTPUT {output}""".format(ref=REF_file_path,
                                                                                                                                             input_indel=self.input()[1].path,
                                                                                                                                             input_snps=self.input()[0].path,
                                                                                                                                             output=self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########14
class Annovar1(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [CombineVariants(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.merged.vcf', '.merged.av'))

    def run(self):
        cmdline = "%s/convert2annovar.pl %s --includeinfo -format vcf4 > %s" % (
            annovar_pro, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


class Annovar2(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Annovar1(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(
            self.input()[0].path.replace('.merged.av', '.merged.anno.%s_multianno.csv' % genome_version))

    def run(self):
        prefix = self.input()[0].path.rpartition('.merged.av')[0]
        cmdline = "%s/table_annovar.pl %s ~/tools/annovar/humandb/ -buildver %s -protocol %s -operation g,r,r,f,f,f,f,f,f -nastring . --remove --otherinfo --csvout --thread %s --outfile %s --argument '-exonicsplicing -splicing 25',,,,,,,," % (
            annovar_pro, prefix + '.merged.av', genome_version, db_names, annovar_thread,
            prefix + '.merged.anno')
        os.system(cmdline)
        record_cmdline(cmdline + '\n\n\n\n' + '{:#^50}'.format('NORMALLY END pipelines'))


class workflow(luigi.Task):
    x = luigi.Parameter()

    def requires(self):
        samples_IDs = str(self.x).split(',')
        for i in samples_IDs:
            yield Annovar2(sampleID=i)


if __name__ == '__main__':
    luigi.run()

    # python -m luigi --module SomaticPipelines_fast_version workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
