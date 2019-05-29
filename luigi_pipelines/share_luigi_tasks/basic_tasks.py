import os
from os.path import dirname

import luigi

from toolkit import run_cmd, valid_path
from .. import config


class QC_trimmomatic(luigi.Task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)
    dry_run = luigi.BoolParameter(default=False)

    def output(self):
        project_name = self.infodict.get("project_name", "")
        odir = self.infodict.get("odir", "")
        odir = config.trim_fmt.format(base=odir,
                                      PN=project_name)
        sample_name = self.infodict.get("SampleID", '')

        if self.PE2:
            ofile_name1 = os.path.join(config.trim_fmt, "{PE_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                PE1_id=sample_name + "_R1")
            ofile_name2 = os.path.join(config.trim_fmt, "{PE_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                PE1_id=sample_name + "_R2")
            return [luigi.LocalTarget(ofile_name1), luigi.LocalTarget(ofile_name2)]
        else:
            ofile_name1 = os.path.join(config.trim_fmt, "{PE_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                SE_id=sample_name)
            return [luigi.LocalTarget(ofile_name1)]

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)  # auto make output dir
        sample_name = self.infodict.get("SampleID", '')
        input1 = self.PE1
        input2 = self.PE2

        if input2:
            cmdline = "java -jar {trimmomatic_jar} PE -threads {thread} {input1} {input2} {ofile1} {base_out}/{PE1_id}.unpaired.fq.gz {ofile2} {base_out}/{PE2_id}.unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                trimmomatic_jar=config.trimmomatic_jar,
                trimmomatic_dir=dirname(config.trimmomatic_jar),
                input1=input1,
                input2=input2,
                PE1_id=sample_name + "_R1",
                PE2_id=sample_name + "_R2",
                ofile1=self.output()[0].path,
                ofile2=self.output()[1].path,
                base_out=dirname(self.output()[0].path),
                thread=config.trimmomatic_thread)
        else:
            cmdline = "java -jar {trimmomatic_jar} SE -threads {thread} {input1} {ofile} ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                trimmomatic_jar=config.trimmomatic_jar,
                trimmomatic_dir=dirname(config.trimmomatic_jar),
                input1=input1,
                SE_id=sample_name,
                ofile=self.output()[0].path,
                base_out=dirname(self.output()[0].path),
                thread=config.trimmomatic_thread)

        run_cmd(cmdline, dry_run=self.dry_run)


class GenerateSam_pair(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        input1 = self.infodict.get("path_R1", "")
        input2 = self.infodict.get("path_R2", "")
        return QC_trimmomatic(PE1=input1,
                              PE2=input2,
                              dry_run=self.dry_run)

    def output(self):
        sample_name = self.infodict.get("SampleID", '')
        project_name = self.infodict.get("project_name", "")
        odir = self.infodict.get("odir", "")

        return luigi.LocalTarget(
            config.output_fmt.format(path=odir,
                                     PN=project_name,
                                     SN=sample_name) + '.sam')

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        sample_name = self.infodict.get("SampleID", '')
        input_file1 = self.input()[0].path
        if len(self.input()) == 1:
            input_file2 = ''
        else:
            input_file2 = self.input()[1].path

        cmdline = "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {o}".format(
            SN=sample_name,
            REF=config.REF_file_path,
            i1=input_file1,
            i2=input_file2,
            o=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class Convertbam(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return GenerateSam_pair(infodict=self.infodict,
                                dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.sam', '.bam'))

    def run(self):
        cmdline = "samtools view -F 0x100 -bSu %s -o %s" % (self.input().path, self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class sorted_bam(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return Convertbam(infodict=self.infodict,
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
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return sorted_bam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        if config.PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "java -Xmx2g -jar {pircard_jar} MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true" % (
                self.input().path,
                self.output().path,
                dirname(self.output().path))
        run_cmd(cmdline, dry_run=self.dry_run)


#########3
class RealignerTargetCreator(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [MarkDuplicate(infodict=self.infodict,
                              dry_run=self.dry_run),
                sorted_bam(infodict=self.infodict,
                           dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.realign.intervals'))

    def run(self):
        if config.PCR_ON:
            input_f = self.input()[1].path
        else:
            input_f = self.input()[0].path
        cmdline = "java -jar {gatk} -T RealignerTargetCreator -nt {thread} -R {REF} -I {input_f} --known {known_gold_vcf} -o {output_f}".format(gatk=config.gatkv36_path,
                                                                                                                                                thread=config.gatk_thread,
                                                                                                                                                REF=config.REF_file_path,
                                                                                                                                                input_f=input_f,
                                                                                                                                                output_f=self.output().path,
                                                                                                                                                known_gold_vcf=config.known_gold_vcf)

        run_cmd(cmdline, dry_run=self.dry_run)


#########4
class IndelRealigner(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [MarkDuplicate(infodict=self.infodict,
                              dry_run=self.dry_run),
                RealignerTargetCreator(infodict=self.infodict,
                                       dry_run=self.dry_run),
                sorted_bam(infodict=self.infodict,
                           dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam', '.realign.bam'))

    def run(self):
        if config.PCR_ON:
            input_f = self.input()[2].path
        else:
            input_f = self.input()[0].path
        cmdline = "java -Xmx5g -jar {gatk} -T IndelRealigner -R {REF} -I {input_f} -targetIntervals {target_Inter} -o {output_f}".format(gatk=config.gatkv36_path,
                                                                                                                                         REF=config.REF_file_path,
                                                                                                                                         input_f=input_f,
                                                                                                                                         target_Inter=self.input()[1].path,
                                                                                                                                         output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run)


#########5
class BaseRecalibrator(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return IndelRealigner(infodict=self.infodict,
                              dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.realign.bam',
                                                           '.recal_data.table'))

    def run(self):
        cmdline = "java -Xmx4g -jar {gatk} -T BaseRecalibrator -nct {gatk_thread} -R {REF} -I {input_f} -knownSites {db_snp} -knownSites {known_gold_vcf} -o {output_f}".format(
            gatk=config.gatkv36_path,
            gatk_thread=config.gatk_thread,
            REF=config.REF_file_path,
            input_f=self.input().path,
            db_snp=config.db_snp,
            known_gold_vcf=config.known_gold_vcf,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


#########6
class PrintReads(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return [IndelRealigner(infodict=self.infodict,
                               dry_run=self.dry_run),
                BaseRecalibrator(infodict=self.infodict,
                                 dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam',
                                                              '.recal_reads.bam'))

    def run(self):
        cmdline = "java -Xmx4g -jar {gatk} -T PrintReads -R {REF} -I {input_f} -BQSR {recal_base} -o {output_f}".format(gatk=config.gatkv36_path,
                                                                                                                        REF=config.REF_file_path,
                                                                                                                        input_f=self.input()[0].path,
                                                                                                                        recal_base=self.input()[1].path,
                                                                                                                        output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run)
