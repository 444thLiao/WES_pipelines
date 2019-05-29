from os.path import join, dirname

import luigi

from main import *
from toolkit import run_cmd, valid_path


class QC_trimmomatic(luigi.Task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)
    dry_run = luigi.BoolParameter(default=False)

    def output(self):
        sample_name = pfn(self.PE1, 'sample_name')
        project_name = pfn(self.PE1, 'project_name')
        return luigi.LocalTarget(join(base_outpath,
                                      "{PN}_result".format(PN=project_name),
                                      "trim_result",
                                      "{SN}_trimed.log".format(SN=sample_name)))

    def run(self):
        valid_path(self.output().path, check_ofile=1)  # auto make output dir

        input1 = self.PE1
        input2 = self.PE2
        suffix = fq_suffix

        if input2 != None:
            cmdline = "java -jar {trimmomatic_jar} PE -threads {thread} {base_in}/{input1}{input_suffix} {base_in}/{input2}{input_suffix} -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                trimmomatic_jar=trimmomatic_jar,
                trimmomatic_dir=dirname(trimmomatic_jar),
                input_suffix=suffix,
                input1=input1,
                input2=input2,
                base_in=base_inpath,
                base_out=dirname(self.output().path),
                output=self.output().path,
                thread=trimmomatic_thread)
        else:
            cmdline = "java -jar {trimmomatic_jar} SE -threads {thread} {base_in}/{input1}{input_suffix} -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                trimmomatic_jar=trimmomatic_jar,
                trimmomatic_dir=os.path.dirname(trimmomatic_jar),
                input_suffix=suffix,
                input1=input1,
                base_in=base_inpath,
                base_out=dirname(self.output().path),
                output=self.output().path,
                thread=trimmomatic_thread)

        run_cmd(cmdline, dry_run=self.dry_run)


class GenerateSam_pair(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        if Pair_data:
            # if pair end sequencing data input.

            if not self_adjust_fn:
                # todo: wrapper into a unify function. instead of place here.
                # if self_adjust_fn set to True, USE pe1_fmt to auto fill name
                input1 = PE1_fmt.format(input=self.sampleID)
                input2 = PE2_fmt.format(input=self.sampleID)
            else:
                # if not set, use glob to auto find.
                input_list = glob.glob(base_inpath + '/*' + self.sampleID + '*')
                if filter_str:
                    input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
                input1 = [_i for _i in input_list if R1_INDICATOR in _i][0]
                input2 = [_i for _i in input_list if R2_INDICATOR in _i][0]
            return QC_trimmomatic(PE1=os.path.basename(input1),
                                  PE2=os.path.basename(input2),
                                  dry_run=self.dry_run)
        else:
            if not self_adjust_fn:
                input1 = SE_fmt.format(input=self.sampleID)
            else:
                input_list = glob.glob(base_inpath + '/*' + self.sampleID + '*')
                if filter_str:
                    input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
                input1 = [_i for _i in input_list if R1_INDICATOR in _i][0]
            return QC_trimmomatic(PE1=os.path.basename(input1),
                                  dry_run=self.dry_run)

    def output(self):
        sample_name = pfn(self.sampleID, 'sample_name')
        project_name = pfn(self.sampleID, 'project_name')

        return luigi.LocalTarget(
            output_fmt.format(path=base_outpath,
                              PN=project_name,
                              SN=sample_name) + '.sam')

    def run(self):
        valid_path(self.output().path, check_ofile=1)

        if Pair_data:
            input_file1 = join(dirname(self.input().path),
                               PE1_fmt.format(input=self.sampleID), )
            input_file2 = join(dirname(self.input().path),
                               PE2_fmt.format(input=self.sampleID), )
        else:
            input_file1 = join(dirname(self.input().path),
                               SE_fmt.format(input=self.sampleID), )
            input_file2 = ''
        cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {o}".format(
            SN=sample_name,
            REF=REF_file_path,
            i1=input_file1,
            i2=input_file2,
            o=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class Convertbam(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return GenerateSam_pair(sampleID=self.sampleID,
                                dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.sam', '.bam'))

    def run(self):
        cmdline = "samtools view -F 0x100 -bSu %s -o %s" % (self.input().path, self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return Convertbam(sampleID=self.sampleID,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.bam', '_sorted.bam'))

    def run(self):
        cmdline = "samtools sort -m {sort_sam_ram} -f -@ {sort_sam_thread} %s %s" % (self.input().path,
                                                                                     self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run)


#########2
class MarkDuplicate(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return sorted_bam(sampleID=self.sampleID,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        if PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "java -Xmx2g -jar {pircard_jar} MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true" % (
                self.input().path,
                self.output().path,
                dirname(self.output().path))
        run_cmd(cmdline, dry_run=self.dry_run)


#########3
class RealignerTargetCreator(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID, dry_run=self.dry_run),
                sorted_bam(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.realign.intervals'))

    def run(self):
        if PCR_ON:
            input_f = self.input()[1].path
        else:
            input_f = self.input()[0].path
        cmdline = "java -jar {gatk} -T RealignerTargetCreator -nt {thread} -R {REF} -I {input_f} --known {known_gold_vcf} -o {output_f}".format(gatk=gatkv36_path,
                                                                                                                                                thread=gatk_thread,
                                                                                                                                                REF=REF_file_path,
                                                                                                                                                input_f=input_f,
                                                                                                                                                output_f=self.output().path,
                                                                                                                                                known_gold_vcf=known_gold_vcf)

        run_cmd(cmdline, dry_run=self.dry_run)


#########4
class IndelRealigner(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID, dry_run=self.dry_run), RealignerTargetCreator(sampleID=self.sampleID, dry_run=self.dry_run),
                sorted_bam(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.realign.bam'))

    def run(self):
        if PCR_ON:
            input_f = self.input()[2].path
        else:
            input_f = self.input()[0].path
        cmdline = "java -Xmx5g -jar {gatk} -T IndelRealigner -R {REF} -I {input_f} -targetIntervals {target_Inter} -o {output_f}".format(gatk=gatkv36_path,
                                                                                                                                         REF=REF_file_path,
                                                                                                                                         input_f=input_f,
                                                                                                                                         target_Inter=self.input()[1].path,
                                                                                                                                         output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run)


#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return IndelRealigner(sampleID=self.sampleID,
                              dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.realign.bam',
                                                              '.recal_data.table'))

    def run(self):
        cmdline = "java -Xmx4g -jar {gatk} -T BaseRecalibrator -nct {gatk_thread} -R {REF} -I {input_f} -knownSites {db_snp} -knownSites {known_gold_vcf} -o {output_f}".format(
            gatk=gatkv36_path,
            gatk_thread=gatk_thread,
            REF=REF_file_path,
            input_f=self.input().path,
            db_snp=db_snp,
            known_gold_vcf=known_gold_vcf,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


#########6
class PrintReads(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID, dry_run=self.dry_run),
                BaseRecalibrator(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam',
                                                              '.recal_reads.bam'))

    def run(self):
        cmdline = "java -Xmx4g -jar {gatk} -T PrintReads -R {REF} -I {input_f} -BQSR {recal_base} -o {output_f}".format(gatk=gatkv36_path,
                                                                                                                        REF=REF_file_path,
                                                                                                                        input_f=self.input()[0].path,
                                                                                                                        recal_base=self.input()[1].path,
                                                                                                                        output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run)
