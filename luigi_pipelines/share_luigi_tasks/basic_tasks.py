import os
from os.path import dirname

import luigi

from api.cal_Cov_script_version import bam2info
from luigi_pipelines import config
from toolkit import run_cmd, valid_path

class base_luigi_task(luigi.Task):
    infodict = luigi.DictParameter(default=dict())
    dry_run = luigi.BoolParameter(default=False)

    def get_log_path(self):
        base_log_path = self.infodict.get("log_path",None)
        if base_log_path is not None:
            return base_log_path



class QC_trimmomatic(base_luigi_task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)

    def output(self):
        odir = self.infodict.get("odir", '')
        project_name = self.infodict.get("project_name", '')
        sample_name = self.infodict.get("SampleID", '')
        if self.PE2:
            ofile_name1 = os.path.join(config.trim_fmt,
                                       "{PE1_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                PE1_id=sample_name + "_R1")
            ofile_name2 = os.path.join(config.trim_fmt,
                                       "{PE1_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                PE1_id=sample_name + "_R2")
            return [luigi.LocalTarget(ofile_name1),
                    luigi.LocalTarget(ofile_name2)]
        else:
            ofile_name1 = os.path.join(config.trim_fmt,
                                       "{PE1_id}.clean.fq.gz").format(
                base=odir,
                PN=project_name,
                PE1_id=sample_name + "_R1")
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
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)

        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())


class GenerateSam_pair(base_luigi_task):

    def requires(self):
        input1 = self.infodict.get("path_R1", "")
        input2 = self.infodict.get("path_R2", "")
        return QC_trimmomatic(PE1=input1,
                              PE2=input2,
                              infodict=self.infodict,
                              dry_run=self.dry_run)

    def output(self):
        odir = self.infodict.get("odir", '')
        project_name = self.infodict.get("project_name", '')
        sample_name = self.infodict.get("SampleID", '')
        return luigi.LocalTarget(config.output_fmt.format(
            path=odir,
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

        cmdline = "bwa mem -M -t {bwa_thread} -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {ofile}".format(
            bwa_thread=config.bwa_thread,
            SN=sample_name,
            REF=config.REF_file_path,
            i1=input_file1,
            i2=input_file2,
            ofile=self.output().path)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path,
                    dry_run=False)

        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())


class Convertbam(base_luigi_task):

    def requires(self):
        return GenerateSam_pair(infodict=self.infodict,
                                dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.sam', '.bam'))

    def run(self):
        if config.samtools_version > 1:
            cmdline = "samtools view -G 100 -bu %s -o %s" % (self.input().path, self.output().path)
        else:
            cmdline = "samtools view -F 0x100 -bSu %s -o %s" % (self.input().path, self.output().path)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


class sorted_bam(base_luigi_task):

    def requires(self):
        return Convertbam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.bam', '_sorted.bam'))

    def run(self):
        if config.samtools_version > 1:
            cmdline = "{samtools_pth} sort -m {sort_sam_ram} -@ {sort_sam_thread} -o {output_path} {input_path} ".format(
                samtools_pth=config.samtools_pro,
                sort_sam_ram=config.sort_sam_ram,
                sort_sam_thread=config.sort_sam_thread,
                output_path=self.output().path,
                input_path=self.input().path)
        else:
            cmdline = "{samtools_pth} sort -m {sort_sam_ram} -f -@ {sort_sam_thread}  {input_path} {output_path} ".format(
                samtools_pth=config.samtools_pro,
                sort_sam_ram=config.sort_sam_ram,
                sort_sam_thread=config.sort_sam_thread,
                output_path=self.output().path,
                input_path=self.input().path)

        run_cmd(cmdline,
                dry_run=self.dry_run, log_file=self.get_log_path())
        cmdline = '{samtools_pth} index {ofile}'.format(samtools_pth=config.samtools_pro,
                                                        ofile=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########2
class MarkDuplicate(base_luigi_task):

    def requires(self):
        return sorted_bam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        valid_path(self.output().path,
                   check_ofile=1)
        if config.PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "java {java_option} -jar {pircard_jar} MarkDuplicates INPUT={input} OUTPUT={output} METRICS_FILE={odir}/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true".format(
                pircard_jar=config.pircard_jar,
                java_option=config.java_option,
                input=self.input().path,
                output=self.output().path,
                odir=os.path.dirname(self.output().path))
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########3
class RealignerTargetCreator(base_luigi_task):

    def requires(self):
        return [MarkDuplicate(infodict=self.infodict,
                              dry_run=self.dry_run),
                sorted_bam(infodict=self.infodict,
                           dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.realign.intervals'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if config.PCR_ON:
            input_f = self.input()[1].path
        else:
            input_f = self.input()[0].path
        cmdline = "java -jar {gatk} -T RealignerTargetCreator -nt {thread} -R {REF} -I {input_f} --known {known_gold_vcf} -o {output_f}".format(
            gatk=config.gatkv36_path,
            thread=config.gatk_thread,
            REF=config.REF_file_path,
            input_f=input_f,
            output_f=self.output().path,
            known_gold_vcf=config.known_gold_vcf)

        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########4
class IndelRealigner(base_luigi_task):

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
        cmdline = "java {java_option} -jar {gatk} -T IndelRealigner -R {REF} -I {input_f} -targetIntervals {target_Inter} -o {output_f}".format(
            gatk=config.gatkv36_path,
            java_option=config.java_option,
            REF=config.REF_file_path,
            input_f=input_f,
            target_Inter=self.input()[1].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        cmdline = 'samtools index %s' % self.output().path
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########5
class BaseRecalibrator(base_luigi_task):

    def requires(self):
        return IndelRealigner(infodict=self.infodict,
                              dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.realign.bam',
                                                           '.recal_data.table'))

    def run(self):
        cmdline = "java {java_option} -jar {gatk} -T BaseRecalibrator -nct {gatk_thread} -R {REF} -I {input_f} -knownSites {db_snp} -knownSites {known_gold_vcf} -o {output_f}".format(
            gatk=config.gatkv36_path,
            java_option=config.java_option,
            gatk_thread=config.gatk_thread,
            REF=config.REF_file_path,
            input_f=self.input().path,
            db_snp=config.db_snp,
            known_gold_vcf=config.known_gold_vcf,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########6
class PrintReads(base_luigi_task):

    def requires(self):
        return [IndelRealigner(infodict=self.infodict,
                               dry_run=self.dry_run,
                               ),
                BaseRecalibrator(infodict=self.infodict,
                                 dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace(".recal_data.table",
                                                              ".recal_reads.bam"))

    def run(self):
        cmdline = "java {java_option} -jar {gatk} -T PrintReads -R {REF} -I {input_f} -BQSR {recal_base} -o {output_f}".format(
            gatk=config.gatkv36_path,
            java_option=config.java_option,
            REF=config.REF_file_path,
            input_f=self.input()[0].path,
            recal_base=self.input()[1].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        cmdline = '{samtools_pth} index {ofile}'.format(samtools_pth=config.samtools_pro,
                                                        ofile=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


############################################################
class cal_coverage_info(base_luigi_task):

    def requires(self):
        return [PrintReads(infodict=self.infodict, dry_run=self.dry_run),
                sorted_bam(infodict=self.infodict, dry_run=self.dry_run)]

    def output(self):
        return [luigi.LocalTarget(self.input()[0].path.replace(".recal_reads.bam",
                                                               ".recal_cov.info")),
                luigi.LocalTarget(self.input()[1].path.replace("_sorted.bam",
                                                               ".sorted_cov.info"))
                ]

    def run(self):
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)

        for _o, _i in zip(self.output(), self.input()):
            if not self.dry_run:
                bam2info(bam_path=_i.path,
                         output_cov=_o.path,
                         bed_file=config.bed_file_path,
                         REF_file=config.REF_file_path
                         )
            else:
                run_cmd("run bam2info for %s" % _i.path, dry_run=self.dry_run, log_file=self.get_log_path())

############################################################
# For quality assessment
class quality_assessment(base_luigi_task):
    """
    calculating the coverage of given bam files.
    Embedded into WES pipelines. Add it into main entry because it is independent.
    """
    tab_file = luigi.Parameter()
    odir = luigi.Parameter()

    def requires(self):
        if "Normal" in self.infodict:
            # if given a pair dict contains both normal and tumor
            return [cal_coverage_info(infodict=self.infodict["Normal"],
                                     dry_run=self.dry_run),
                    cal_coverage_info(infodict=self.infodict["Tumor"],
                                      dry_run=self.dry_run)]
        else:
            # just a single sample dict.
            return cal_coverage_info(infodict=self.infodict,
                                 dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              'quality_accessment_raw.csv'))

    def run(self):
        # summarize python script
        # it will iterate all samples contains at `self.tab_file`
        py_file = os.path.join(dirname(dirname(dirname(__file__))),
                               "api",
                               "quality_accessment.py")

        run_cmd("python3 {pyfile} -i {input} -o {output}".format(
            pyfile=py_file,
            input=self.tab_file,
            output=self.odir),
            dry_run=self.dry_run,
            log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)
