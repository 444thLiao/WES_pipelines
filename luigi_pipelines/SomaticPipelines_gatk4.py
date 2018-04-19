import glob, time
import luigi
from parse_file_name import pfn
from collections import defaultdict
from main import *


def record_cmdline(message, default=base_outpath + '/somatic_pipelines.log'):
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
        project_name = pfn(self.PE1, 'project_name')
        output1 = PE1_fmt.format(input=pfn(self.PE1, 'sample_name'))
        return luigi.LocalTarget(
            '{base}/{PN}_result/trim_result/{output1}.clean.fq.gz'.format(base=base_outpath, PN=project_name,
                                                                          output1=output1))

    def run(self):
        sample_name = pfn(self.PE1, 'sample_name')
        project_name = pfn(self.PE1, 'project_name')
        log_name = '{base}/{PN}_result/trim_result/{SN}_trimed.log'.format(base=base_outpath, PN=project_name,
                                                                           SN=sample_name)
        if not os.path.isdir('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name)):
            os.makedirs('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name))

        input1 = self.PE1
        input2 = self.PE2
        output1 = PE1_fmt.format(input=pfn(self.PE1, 'sample_name'))
        output2 = PE2_fmt.format(input=pfn(self.PE2, 'sample_name'))

        if input2:
            cmdline = "java -jar /home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{output1}.clean.fq.gz {base_out}/{output1}.unpaired.fq.gz {base_out}/{output2}.clean.fq.gz {base_out}/{output2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                input1=input1, input2=input2, base_in=base_inpath, base_out=os.path.dirname(log_name),
                output1=output1, output2=output2,
                output=log_name)
            os.system(cmdline)
            record_cmdline(cmdline)
        else:
            cmdline = "java -jar /home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 20 {base_in}/{input1}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                input1=input1, base_in=base_inpath, base_out=os.path.dirname(log_name),
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
        return luigi.LocalTarget(self.input()[0].path.replace('.sam', '.bam'))

    def run(self):
        cmdline = "%s view -@ 20 -F 0x100 -T %s -b %s -o %s" % (
        samtools_pro, REF_file_path, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)
        cmdline = "touch %s" % self.input()[0].path
        os.system(cmdline)
        record_cmdline(cmdline)


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Convertbam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.bam', '_sorted.bam'))

    def run(self):
        cmdline = "%s sort -m 2G -@ 30 %s -o %s" % (samtools_pro, self.input()[0].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)
        cmdline = '%s index -@ 30 %s' % (samtools_pro, self.output().path)
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
            cmdline = "%s MarkDuplicates --java-options '-Xmx30g' --INPUT %s --OUTPUT %s --METRICS_FILE %s/dedup_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true -AS true" % (
                gatk_pro, self.input()[0].path, self.output().path, self.output().path.rpartition('/')[0])
        os.system(cmdline)
        record_cmdline(cmdline)


#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.recal_data.table'))

    def run(self):
        cmdline = "%s BaseRecalibrator --java-options '-Xmx30g' --reference %s --input %s --known-sites %s --known-sites %s --output %s" % (
            gatk_pro, REF_file_path, self.input()[0].path, db_snp, known_gold_cvf, self.output().path)
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
        cmdline = "%s ApplyBQSR --java-options '-Xmx30g' --reference %s --input %s --bqsr-recal-file %s --output %s" % (
            gatk_pro, REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)
        cmdline = '%s index -@ 30 %s' % (samtools_pro, self.output().path)
        os.system(cmdline)
        record_cmdline(cmdline)


#########somatic pipeline
class MuTect2_pair(luigi.Task):
    ####combine N & T
    sample_IDs = luigi.Parameter()

    def requires(self):
        sampleIDs = self.sample_IDs.split(',')
        # TODO
        # process situation which is 'len(sampleIDs) >2'. It may be some sample need to combine 2 or 3 or more.
        return [PrintReads(sampleID=sampleIDs[0]), PrintReads(sampleID=sampleIDs[1])]

    def output(self):
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = [k for k, v in pair_bucket.items() if set(v) == set(sampleIDs)][0]

        output_path = somatic_pair_output_fmt.format(path=base_outpath, PN=Project_ID, PairN=pair_name) + '.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        sampleIDs = self.sample_IDs.split(',')

        output_dir = self.output().path.rpartition('/')[0]

        if os.path.isdir(output_dir) != True:
            os.makedirs(output_dir)
        input_tumor = ''
        input_normal = ''
        normal_name = ''
        tumor_name = ''

        if pfn(sampleIDs[0], 'mt2_for') == NORMAL_SIG:
            input_normal = self.input()[0].path
            input_tumor = self.input()[1].path
            normal_name = pfn(sampleIDs[0], 'sample_name')
            tumor_name = pfn(sampleIDs[1], 'sample_name')
        elif pfn(sampleIDs[0], 'mt2_for') == TUMOR_SIG:
            input_normal = self.input()[1].path
            input_tumor = self.input()[0].path
            normal_name = pfn(sampleIDs[1], 'sample_name')
            tumor_name = pfn(sampleIDs[0], 'sample_name')

        prefix = self.output().path.rpartition('.bam')[0]
        if bed_file_path:
            suffix_str = " --intervals %s" % bed_file_path
        else:
            suffix_str = ''
        cmdline = "{gatk} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_normal} -normal {N_name} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam".format(
            REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_tumor=input_tumor, input_normal=input_normal,
            gatk=gatk_pro,
            N_name=normal_name, T_name=tumor_name, prefix=prefix) + suffix_str
        os.system(cmdline)
        record_cmdline(cmdline)


class MuTect2_single(luigi.Task):
    sample_NT = luigi.Parameter()

    def requires(self):
        return [PrintReads(sampleID=self.sample_NT)]

    def output(self):

        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')

        output_path = somatic_single_output_fmt.format(path=base_outpath, PN=Project_ID, SN=sample_name) + '.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        input1 = self.input()[0].path
        mt2_id = pfn(self.sample_NT, 'mt2_for')
        prefix = self.output().path.rpartition('.bam')[0]
        output_dir = self.output().path.rpartition('/')[0]

        if os.path.isdir(output_dir) != True:
            os.makedirs(output_dir)

        if bed_file_path:
            suffix_str = " --intervals %s" % bed_file_path
        else:
            suffix_str = ''

        if mt2_id == NORMAL_SIG:
            cmdline = "{gatk} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam".format(
                gatk=gatk_pro,
                REF=REF_file_path, db_snp=db_snp, input_tumor=input1, prefix=prefix,
                T_name=pfn(self.sample_NT, 'sample_name')) + suffix_str
            os.system(cmdline)
            record_cmdline(cmdline)

        # Normal only
        else:
            cmdline = "{gatk} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam --tumor-lod-to-emit 4".format(
                gatk=gatk_pro,
                REF=REF_file_path, db_snp=db_snp, input_tumor=input1, prefix=prefix,
                T_name=pfn(self.sample_NT, 'sample_name')) + suffix_str
            os.system(cmdline)
            record_cmdline(cmdline)
            # Tumor only


class Annovar1(luigi.Task):
    sample_ID = luigi.Parameter()

    def requires(self):
        if self.sample_ID in pair_bucket:
            pair_value = pair_bucket[self.sample_ID]
            return [MuTect2_pair(sample_IDs=','.join(pair_value))]
        else:
            return [MuTect2_single(sample_NT=self.sample_ID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.bam')[0] + '.merged.av')

    def run(self):
        for i in self.input():
            prefix = i.path.rpartition('.bam')[0]
            cmdline = "%s/convert2annovar.pl %s --includeinfo -format vcf4 > %s" % (
                annovar_pro, prefix + '.vcf', prefix + '.merged.av')
            os.system(cmdline)
            record_cmdline(cmdline)


class Annovar2(luigi.Task):
    sample_ID = luigi.Parameter()

    def requires(self):
        return [Annovar1(sample_ID=self.sample_ID)]

    def output(self):
        return luigi.LocalTarget(
            self.input()[0].path.replace('.merged.av', '.merged.anno.%s_multianno.csv' % genome_version))

    def run(self):
        prefix = self.input()[0].path.rpartition('.merged.av')[0]
        cmdline = "%s/table_annovar.pl %s /home/liaoth/tools/annovar/humandb/ -buildver %s -protocol %s -operation g,r,r,f,f,f,f,f,f -nastring . --remove --otherinfo --csvout --thread %s --outfile %s --argument '-exonicsplicing -splicing 25',,,,,,,, " % (
            annovar_pro, prefix + '.merged.av', genome_version, db_names, annovar_thread,
            prefix + '.merged.anno')
        os.system(cmdline)
        record_cmdline(cmdline + '\n\n\n\n' + '{:#^50}'.format('NORMALLY END pipelines'))


class workflow(luigi.Task):
    x = luigi.Parameter()

    def requires(self):
        samples_IDs = str(self.x).split(',')

        pair_bucket = defaultdict(list)
        for _x in samples_IDs:
            pair_bucket[pfn(_x, 'pair_name')].append(_x)
        adjust_multiple = []
        for each in pair_bucket.keys():
            if len(pair_bucket[each]) > 2:
                tmp = pair_bucket[each]
                only_normal = [_ for _ in tmp if pfn(_, 'mt2_for') == NORMAL_SIG][0]
                for _each in tmp:
                    if pfn(_each, 'mt2_for') == TUMOR_SIG and pfn(_each, 'sample_name').replace(TUMOR_SIG, '') != each:
                        adjust_multiple.append((pfn(_each, 'sample_name').replace(TUMOR_SIG, ''), [only_normal, _each]))
                    elif pfn(_each, 'mt2_for') == TUMOR_SIG and pfn(_each, 'sample_name').replace(TUMOR_SIG,
                                                                                                  '') == each:
                        adjust_multiple.append((each, [only_normal, _each]))
        pair_bucket.update(dict(adjust_multiple))
        global pair_bucket
        ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}

        samples_IDs += [_x for _x in pair_bucket.keys()]

        if debug_:
            import pdb;
            pdb.set_trace()
        for i in samples_IDs:
            yield Annovar2(sample_ID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
