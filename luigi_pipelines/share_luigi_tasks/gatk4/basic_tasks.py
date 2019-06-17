from os.path import dirname

import luigi

from toolkit import run_cmd, valid_path
from luigi_pipelines import config
from luigi_pipelines.share_luigi_tasks import sorted_bam

class base_luigi_task(luigi.Task):
    infodict = luigi.DictParameter(default=dict())
    dry_run = luigi.BoolParameter(default=False)

    def get_log_path(self):
        base_log_path = self.infodict.get("log_path",None)
        if base_log_path is not None:
            return base_log_path

#########2
class MarkDuplicate(base_luigi_task):

    def requires(self):
        return sorted_bam(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('_sorted.bam', '.dedup.bam'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if config.PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "{gatk4} MarkDuplicates --java-options '{java_option}' --INPUT {input_f} --OUTPUT {output_f} --METRICS_FILE {odir}/dedup_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true -AS true".format(
                gatk4=config.gatk_pro,
                java_option=config.java_option,
                input_f=self.input().path,
                output_f=self.output().path,
                odir=dirname(self.output().path)
            )

        run_cmd(cmdline, dry_run=self.dry_run,log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


#########5
class BaseRecalibrator(base_luigi_task):

    def requires(self):
        return MarkDuplicate(infodict=self.infodict,
                             dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.dedup.bam',
                                                           '.recal_data.table'))

    def run(self):
        cmdline = "{gatk4} BaseRecalibrator --java-options '{java_option}' --reference {REF} --input {input_f} --known-sites {db_snp} --known-sites {known_gold_vcf} --output {output_f}".format(
            gatk4=config.gatk_pro,
            java_option=config.java_option,
            REF=config.REF_file_path,
            input_f=self.input().path,
            db_snp=config.db_snp,
            known_gold_vcf=config.known_gold_vcf,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run,log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

#########6
class PrintReads(base_luigi_task):

    def requires(self):
        return [MarkDuplicate(infodict=self.infodict, dry_run=self.dry_run),
                BaseRecalibrator(infodict=self.infodict, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.recal_reads.bam'))

    def run(self):
        cmdline = "{gatk4} ApplyBQSR --java-options '{java_option}' --reference {REF} --input {input_f} --bqsr-recal-file {recal_base} --output {output_f}".format(
            gatk4=config.gatk_pro,
            java_option=config.java_option,
            REF=config.REF_file_path,
            input_f=self.input()[0].path,
            recal_base=self.input()[1].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run,log_file=self.get_log_path())
        cmdline = 'samtools index -@ %s %s' % (config.sort_sam_thread,
                                               self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run,log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)