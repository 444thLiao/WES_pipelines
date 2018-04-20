###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample Output a gemini db.
### @GPZ-bioinfo, 20170313
###############################################################################################
import luigi,time
import os
from special_fun import Add_cov_ino_in_vcf as P_vcf

def record_cmdline(message, default=base_outpath + '/%s_pipelines.log' % os.path.basename(__file__.replace('.py',''))):
    if os.path.isfile(default):
        with open(default, 'a') as f1:
            f1.write(time.ctime() + ' ' * 4 + message + '\n')
    else:
        with open(default, 'w') as f1:
            f1.write('{:#^40}'.format('Starting the %s pipelines.' % os.path.basename(__file__.replace('.py',''))))
            f1.write(time.ctime() + ' ' * 4 + message + '\n')


class QC_trimmomatic(luigi.Task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)

    def output(self):
        project_name = pfn(self.PE1, 'project_name')
        output1 = PE1_fmt.format(input=pfn(self.PE1, 'sample_name'))
        return luigi.LocalTarget(os.path.join(trim_fmt.format(base=base_outpath, PN=project_name),
                                              '/%s.clean.fq.gz' % output1))

    def run(self):
        sample_name = pfn(self.PE1, 'sample_name')
        project_name = pfn(self.PE1, 'project_name')
        trim_r_path = trim_fmt.format(base=base_outpath, PN=project_name)
        log_name = os.path.join(trim_r_path, '%s_trimed.log' % sample_name)

        if not os.path.isdir(trim_r_path):
            os.makedirs(trim_r_path)

        input1 = self.PE1
        input2 = self.PE2
        output1 = PE1_fmt.format(input=pfn(self.PE1, 'sample_name'))
        output2 = PE2_fmt.format(input=pfn(self.PE2, 'sample_name'))

        if input2:
            cmdline = "java -jar {trimmomatic_jar} PE -threads 20 {base_in}/{input1}{fq_suffix} {base_in}/{input2}{fq_suffix} -trimlog {output} {base_out}/{output1}.clean.fq.gz {base_out}/{output1}.unpaired.fq.gz {base_out}/{output2}.clean.fq.gz {base_out}/{output2}.unpaired.fq.gz ILLUMINACLIP:{trimmomatic_jar_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                trimmomatic_jar=trimmomatic_jar, trimmomatic_jar_dir=os.path.dirname(trimmomatic_jar),
                input1=input1, input2=input2, base_in=base_inpath, base_out=os.path.dirname(log_name),
                output1=output1, output2=output2, fq_suffix=fq_suffix,
                output=log_name)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            cmdline = "java -jar {trimmomatic_jar} SE -threads 20 {base_in}/{input1}{fq_suffix} -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:{trimmomatic_jar_dir}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                trimmomatic_jar=trimmomatic_jar, trimmomatic_jar_dir=os.path.dirname(trimmomatic_jar),
                input1=input1, base_in=base_inpath, base_out=os.path.dirname(log_name), fq_suffix=fq_suffix,
                output=log_name)
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
                input1 = [_i.replace(fq_suffix, '') for _i in input_list if R1_INDICATOR in _i][0]
                input2 = [_i.replace(fq_suffix, '') for _i in input_list if R2_INDICATOR in _i][0]
            return QC_trimmomatic(PE1=os.path.basename(input1), PE2=os.path.basename(input2))
        else:
            if not self_adjust_fn:
                input1 = SE_fmt.format(input=self.sampleID)
            else:
                input_list = glob.glob(base_inpath + '/*' + self.sampleID + '*')
                if filter_str:
                    input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
                input1 = [_i.replace(fq_suffix, '') for _i in input_list if R1_INDICATOR in _i][0]
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
            input1 = self.input().path
            input2 = self.input().path.replace(R1_INDICATOR, R2_INDICATOR)
            if not os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)):
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2} > {o}".format(
                SN=sample_name, REF=REF_file_path, i1=input1, i2=input2, o=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            input1 = self.input().path
            if not os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)):
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} > {o}".format(
                SN=sample_name, REF=REF_file_path, i1=input1, o=self.output().path)
            os.system(cmdline)
            record_cmdline(cmdline)

class Convertbam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [GenerateSam_pair(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.sam','.bam'))

    def run(self):
        os.system("samtools view -F 0x100 -bSu %s -o %s" % (self.input()[0].path, self.output().path))


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Convertbam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.bam','_sorted.bam'))

    def run(self):
        os.system("samtools sort -m 100G -f -@ 30 %s %s" % (self.input()[0].path, self.output().path))
        os.system('samtools index %s' % self.output().path)

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
            cmdline = "java -Xmx2g -jar ~/tools/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true" % (
            self.input()[0].path, self.output().path, self.output().path.rpartition('/')[0])
        os.system(cmdline)

#########3
class RealignerTargetCreator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID), sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.realign.intervals'))

    def run(self):
        if PCR_ON:
            cmdline = "java -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R %s -I %s --known %s -o %s" % (
                REF_file_path, self.input()[1].path, known_gold_cvf, self.output().path)

        else:
            cmdline = "java -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R %s -I %s --known %s -o %s" % (
                REF_file_path, self.input()[0].path, known_gold_cvf, self.output().path)
        os.system(cmdline)

#########4
class IndelRealigner(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID), RealignerTargetCreator(sampleID=self.sampleID),
                sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.realign.bam'))

    def run(self):
        if PCR_ON:
            cmdline = "java -Xmx5g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (
                REF_file_path, self.input()[2].path, self.input()[1].path, self.output().path)
        else:
            cmdline = "java -Xmx5g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path)
        os.system(cmdline)
        cmdline = 'samtools index %s' % self.output().path
        os.system(cmdline)

#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam','.recal_data.table'))

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
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam','.recal_reads.bam'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads -R %s -I %s -BQSR %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


#########7
class HaplotypeCaller(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [PrintReads(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.recal_reads.bam','.raw_variants.vcf'))

    def run(self):
        if bed_file_path != '':

            os.system(
                "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 30 -R %s -I %s -L %s --genotyping_mode DISCOVERY --dbsnp %s -stand_call_conf 10 -stand_emit_conf 5 -A AlleleBalance -A Coverage -A FisherStrand -o %s" % (
                    REF_file_path, self.input()[0].path, bed_file_path, db_snp, self.output().path))
        else:
            os.system(
                "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 30 -R %s -I %s --genotyping_mode DISCOVERY --dbsnp %s -stand_call_conf 10 -stand_emit_conf 5 -A AlleleBalance -A Coverage -A FisherStrand -o %s" % (
                    REF_file_path, self.input()[0].path, db_snp, self.output().path))


#########9
class SelectVariants_a(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [HaplotypeCaller(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_variants.vcf','.raw_snps.vcf'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R %s -V %s -selectType SNP -o %s" % (
                REF_file_path, self.input()[0].path, self.output().path))


#########10
class VariantFiltration_a(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [SelectVariants_a(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_snps.vcf','.filter_snps.vcf'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\" -o %s" % (
                REF_file_path, self.input()[0].path, self.output().path))


#########11
class SelectVariants_b(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [HaplotypeCaller(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_variants.vcf','.raw_indels.vcf'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R %s -V %s -selectType INDEL -o %s" % (
                REF_file_path, self.input()[0].path, self.output().path))


#########12
class VariantFiltration_b(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [SelectVariants_b(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_indels.vcf','.filter_indels.vcf'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o %s" % (
                REF_file_path, self.input()[0].path, self.output().path))


#########13
class CombineVariants(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [VariantFiltration_a(sampleID=self.sampleID), VariantFiltration_b(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace('.filter_indels.vcf','.merged.vcf'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T CombineVariants -R %s --variant:indel %s --variant:snp %s --interval_padding 25 --out %s --setKey set --genotypemergeoption UNSORTED" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))


#########14
class Add_cov_infos(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [CombineVariants(sampleID=self.sampleID),
                HaplotypeCaller(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.merged.vcf', '.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_SO(self.input()[1].path, self.input()[0].path, self.output().path, REF_file_path)


#########15
class vt_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Add_cov_infos(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.added_cov.vcf', '.vt.vcf'))

    def run(self):
        os.system("{vt} decompose -s {input_vcf} | {vt} normalize -r {REF}- > {vt_vcf}".format(
            vt=vt_pro, input_vcf=self.input()[0].path, REF=REF_file_path, vt_vcf=self.output().path
        ))


class vep_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [vt_part(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.vt.vcf', '.vep.vcf.gz'))

    def run(self):
        os.system("""{vep} -i {vt_vcf} --cache --merged --fasta {REF} --sift b --polyphen b --symbol --numbers --biotype
        --total_length --canonical --ccds -o {vep_output_vcf} --vcf --hgvs --gene_phenotype --uniprot
        --force_overwrite --port 3337 --domains --regulatory --protein --tsl --variant_class --fork {threads} --force
        --no_stats >> {vep_log} 2>&1""".format(vep=vep_pro,
                                               vt_vcf=self.input()[0].path,
                                               REF=REF_file_path,
                                               vep_output_vcf=self.input()[0].path.replace('.vt.vcf', '.vep.vcf'),
                                               threads=20,
                                               vep_log=self.input()[0].path.replace('.vt.vcf', '.vep.log')))
        os.system('{bgzip} -c {vep_output_vcf} > {vep_output_vcf_gz}'.format(bgzip=bgzip_pro,
                                                                             vep_output_vcf=self.input()[
                                                                                 0].path.replace('.vt.vcf', '.vep.vcf'),
                                                                             vep_output_vcf_gz=self.output().path))
        os.system('{tabix} -p vcf {vep_output_vcf_gz}'.format(tabix=tabix_pro,
                                                              vep_output_vcf_gz=self.output().path))


class gemini_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [vep_part(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.vep.vcf.gz', '.gemini.db'))

    def run(self):
        os.system("""{gemini} load --cores {threads} -t VEP -v {vep_output_vcf_gz} {Output_db};
        {gemini} annotate -f {vep_output_vcf_gz} -a extract
        -c SAD,SAF,AF,AD,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR
        -t text,float,float,text,float,float,float,float,float
        -o list,list,list,list,mean,mean,mean,mean,mean {Output_db} >> {gemini_log} 2>&1""".format(
            gemini=gemini_pro,
            threads=20,
            vep_output_vcf_gz=self.input()[0].path,
            Output_db=self.output().path,
            gemini_log=self.output().path + '.log'))



class workflow(luigi.Task):
    x = luigi.Parameter()

    def requires(self):
        samples_IDs = str(self.x).split(',')
        for i in samples_IDs:
            yield gemini_part(sampleID=i)


if __name__ == '__main__':
    luigi.run()

from main import *
# python -m luigi --module SomaticPipelines_fast_version workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
