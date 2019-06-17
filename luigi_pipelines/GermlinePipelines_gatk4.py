###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample 
### @GPZ-bioinfo, 20190525
###############################################################################################

import luigi

from luigi_pipelines import config, run_cmd, valid_path
from luigi_pipelines.GermlinePipelines import HaplotypeCaller, SelectVariants, VariantFiltration, CombineVariants
from luigi_pipelines.share_luigi_tasks import Annovar1, Annovar2
from luigi_pipelines.share_luigi_tasks.gatk4 import PrintReads


#########7

class HaplotypeCaller(HaplotypeCaller):
    def requires(self):
        return PrintReads(infodict=self.infodict, dry_run=self.dry_run)

    def run(self):
        valid_path(self.output().path, check_ofile=1)

        if config.bed_file_path != '':
            extra_str = " --intervals {}".format(config.bed_file_path)
        else:
            extra_str = ""
        cmdline = "{gatk4} HaplotypeCaller --java-options '-Xmx30g' --native-pair-hmm-threads 30 --reference {ref} --input {input} --genotyping-mode DISCOVERY --dbsnp {dbsnp} -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {output} {extra_str}".format(
            ref=config.REF_file_path,
            input=self.input().path,
            dbsnp=config.db_snp,
            output=self.output().path,
            extra_str=extra_str,
            gatk4=config.gatk_pro)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


#########9
class SelectVariants(SelectVariants):

    def requires(self):
        return HaplotypeCaller(infodict=self.infodict,
                               dry_run=self.dry_run)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            selecttype = "SNP"
        elif self.object_type == "indel":
            selecttype = "INDEL"
        else:
            raise Exception

        cmdline = "{gatk4} SelectVariants --java-options '-Xmx4g' -R {REF} -V {input_f} -select-type {selecttype} -O {output_f}".format(
            gatk4=config.gatk_pro,
            REF=config.REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            selecttype=selecttype)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


#########10
class VariantFiltration(VariantFiltration):
    def requires(self):
        return SelectVariants(infodict=self.infodict,
                              dry_run=self.dry_run,
                              object_type=self.object_type)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            filterExpression = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        elif self.object_type == "indel":
            filterExpression = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        else:
            raise Exception

        cmdline = """{gatk4} VariantFiltration --java-options '-Xmx4g' -R {REF} -V {input_f} --filter-expression "{filterExpression}" --filter-name \"my_{object_type}_filter\" -O {output_f}""".format(
            gatk4=config.gatk_pro,
            REF=config.REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            filterExpression=filterExpression,
            object_type=self.object_type)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


#########13
class CombineVariants(CombineVariants):

    def requires(self):
        required_task = {ot: VariantFiltration(infodict=self.infodict,
                                               dry_run=self.dry_run,
                                               object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        cmdline = """{gatk4} MergeVcfs --java-options "-Xmx4g" -R {REF} --INPUT {input_indel} --INPUT {input_snp} --OUTPUT {output_f}""".format(
            gatk4=config.gatk_pro,
            REF=config.REF_file_path,
            input_indel=self.input()["indel"].path,
            input_snp=self.input()["snp"].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())


class new_Annovar1(Annovar1):
    def requires(self):
        return CombineVariants(infodict=self.infodict,
                               dry_run=self.dry_run)


class new_Annovar2(Annovar2):
    def requires(self):
        return [new_Annovar1(infodict=self.infodict,
                             dry_run=self.dry_run)]


if __name__ == '__main__':
    luigi.run()

    #
