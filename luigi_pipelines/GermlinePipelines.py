###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample 
### @GPZ-bioinfo, 20190525
###############################################################################################

from os.path import dirname

import luigi
import os

from luigi_pipelines import config, run_cmd, valid_path
from luigi_pipelines.share_luigi_tasks import PrintReads, Annovar1, Annovar2




#########7
class HaplotypeCaller(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return PrintReads(infodict=self.infodict, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.recal_reads.bam',
                                                           '.raw_variants.vcf'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)

        if config.bed_file_path != '':
            extra_str = "-L %s" % config.bed_file_path
        else:
            extra_str = ''
        cmdline = "java -Xmx4g -jar {gatk} -T HaplotypeCaller -nct {gatk_thread} -R {REF} -I {input} {extra_str} --genotyping_mode DISCOVERY --dbsnp {db_snp} -stand_call_conf 10 -stand_emit_conf 5 -A AlleleBalance -A Coverage -A FisherStrand -o {output_f}".format(
            gatk=config.gatkv36_path,
            gatk_thread=config.gatk_thread,
            REF=config.REF_file_path,
            input=self.input().path,
            extra_str=extra_str,
            db_snp=config.db_snp,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


#########9
class SelectVariants(luigi.Task):
    infodict = luigi.DictParameter()
    object_type = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return HaplotypeCaller(infodict=self.infodict,
                               dry_run=self.dry_run)

    def output(self):
        if self.object_type == "snp":
            ofile_name = '.raw_snps.vcf'
        elif self.object_type == "indel":
            ofile_name = '.raw_indels.vcf'
        else:
            raise Exception

        return luigi.LocalTarget(self.input().path.replace('.raw_variants.vcf',
                                                           ofile_name))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            selecttype = "SNP"
        elif self.object_type == "indel":
            selecttype = "INDEL"
        else:
            raise Exception

        cmdline = "java -Xmx4g -jar {gatk} -T SelectVariants -R {REF} -V {input_f} -selectType {selecttype} -o {output_f}".format(
            gatk=config.gatkv36_path,
            REF=config.REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            selecttype=selecttype)
        run_cmd(cmdline, dry_run=self.dry_run)


#########10
class VariantFiltration(luigi.Task):
    infodict = luigi.DictParameter()
    object_type = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return SelectVariants(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.raw_',
                                                           '.filter_'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            filterExpression = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        elif self.object_type == "indel":
            filterExpression = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        else:
            raise Exception
        cmdline = """java -Xmx4g -jar {gatk} -T VariantFiltration -R {REF} -V {input_f} --filterExpression "{filterExpression}" --filterName \"my_{object_type}_filter\" -o {output_f}""".format(
            gatk=config.gatkv36_path,
            REF=config.REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            filterExpression=filterExpression,
            object_type=self.object_type)
        run_cmd(cmdline, dry_run=self.dry_run)


#########13
class CombineVariants(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        required_task = {ot: VariantFiltration(infodict=self.infodict,
                                               dry_run=self.dry_run,
                                               object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task

    def output(self):
        return luigi.LocalTarget(self.input()["indel"].path.replace('.filter_indels.vcf',
                                                                    '.merged.vcf'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        cmdline = "java -Xmx4g -jar {gatk} -T CombineVariants -R {REF} --variant:indel {input_indel} --variant:snp {input_snp} --interval_padding 25 --out {output_f} --setKey set --genotypemergeoption UNSORTED".format(
            gatk=config.gatkv36_path,
            REF=config.REF_file_path,
            input_indel=self.input()["indel"].path,
            input_snp=self.input()["snp"].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class new_Annovar1(Annovar1):
    def requires(self):
        return CombineVariants(infodict=self.infodict,
                               dry_run=self.dry_run)


class new_Annovar2(Annovar2):
    def requires(self):
        return new_Annovar1(infodict=self.infodict,
                            dry_run=self.dry_run)


if __name__ == '__main__':
    luigi.run()

    # python3 GermlinePipelines.py new_Annovar2
