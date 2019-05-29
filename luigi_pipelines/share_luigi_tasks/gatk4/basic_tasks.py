from os.path import dirname

import luigi

from main import *
from toolkit import run_cmd, valid_path
from .. import sorted_bam


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
        valid_path(self.output().path, check_ofile=1)
        if PCR_ON:
            cmdline = "touch %s" % self.output().path
        else:
            cmdline = "{gatk4} MarkDuplicates --java-options '-Xmx30g' --INPUT {input_f} --OUTPUT {output_f} --METRICS_FILE {odir}/dedup_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true -AS true".format(
                gatk4=gatk_pro,
                input_f=self.input().path,
                output_f=self.output().path,
                odir=dirname(self.output().path)
            )

        run_cmd(cmdline, dry_run=self.dry_run)


#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return MarkDuplicate(sampleID=self.sampleID,
                             dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.dedup.bam',
                                                           '.recal_data.table'))

    def run(self):
        cmdline = "{gatk4} BaseRecalibrator --java-options '-Xmx30g' --reference {REF} --input {input_f} --known-sites {db_snp} --known-sites {known_gold_vcf} --output {output_f}".format(
            gatk4=gatk_pro,
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
        return [MarkDuplicate(sampleID=self.sampleID, dry_run=self.dry_run),
                BaseRecalibrator(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam',
                                                              '.recal_reads.bam'))

    def run(self):
        cmdline = "{gatk4} ApplyBQSR --java-options '-Xmx30g' --reference {REF} --input {input_f} --bqsr-recal-file {recal_base} --output {output_f}".format(
            gatk4=gatk_pro,
                                                                                                                        REF=REF_file_path,
                                                                                                                        input_f=self.input()[0].path,
                                                                                                                        recal_base=self.input()[1].path,
                                                                                                                        output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)
        cmdline = 'samtools index -@ %s %s' % (sort_sam_thread,
                                               self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)