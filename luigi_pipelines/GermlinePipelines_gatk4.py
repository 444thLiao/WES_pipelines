###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample 
### @GPZ-bioinfo, 20190525
###############################################################################################

import luigi

from luigi_pipelines.share_luigi_tasks import Annovar1, Annovar2
from luigi_pipelines.share_luigi_tasks.gatk4 import PrintReads
from main import *
from .. import valid_path, run_cmd


#########7
class HaplotypeCaller(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return PrintReads(sampleID=self.sampleID, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.recal_reads.bam',
                                                              '.raw_variants.vcf'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)

        if bed_file_path != '':
            extra_str = " --intervals {}".format(bed_file_path)
        else:
            extra_str = ""
        cmdline = "{gatk4} HaplotypeCaller --java-options '-Xmx30g' --native-pair-hmm-threads 30 --reference {ref} --input {input} --genotyping-mode DISCOVERY --dbsnp {dbsnp} -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {output} {extra_str}".format(
            ref=REF_file_path,
            input=self.input()[0].path,
            dbsnp=db_snp,
            output=self.output().path,
            extra_str=extra_str,
            gatk4=gatk_pro)
        run_cmd(cmdline, dry_run=self.dry_run)


#########9
class SelectVariants(luigi.Task):
    sampleID = luigi.Parameter()
    object_type = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return HaplotypeCaller(sampleID=self.sampleID,
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

        cmdline = "{gatk4} SelectVariants --java-options '-Xmx4g' -R {REF} -V {input_f} -select-type {selecttype} -O {output_f}".format(
            gatk4=gatk_pro,
            REF=REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            selecttype=selecttype)
        run_cmd(cmdline, dry_run=self.dry_run)


#########10
class VariantFiltration(luigi.Task):
    sampleID = luigi.Parameter()
    object_type = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return SelectVariants(sampleID=self.sampleID,
                              dry_run=self.dry_run,
                              object_type=self.object_type)

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.raw_',
                                                              '.filter_'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if self.object_type == "snp":
            filterExpression = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        elif self.object_type == "indel":
            filterExpression = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        else:
            raise Exception

        cmdline = """{gatk4} VariantFiltration --java-options '-Xmx4g' -R {REF} -V {input_f} --filter-expression "{filterExpression}" --filter-name \"my_{object_type}_filter\" -O {output_f}""".format(
            gatk4=gatk_pro,
            REF=REF_file_path,
            input_f=self.input().path,
            output_f=self.output().path,
            filterExpression=filterExpression,
            object_type=self.object_type)
        run_cmd(cmdline, dry_run=self.dry_run)


#########13
class CombineVariants(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        required_task = {ot: VariantFiltration(sampleID=self.sampleID,
                                               dry_run=self.dry_run,
                                               object_type=ot)
                         for ot in ["snp", "indel"]}
        return required_task

    def output(self):
        return luigi.LocalTarget(self.input()["indel"].path.replace('.filter_indels.vcf',
                                                                    '.merged.vcf'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        cmdline = """{gatk4} MergeVcfs --java-options "-Xmx4g" -R {REF} --INPUT {input_indel} --INPUT {input_snp} --OUTPUT {output_f}""".format(
            gatk4=gatk_pro,
            REF=REF_file_path,
            input_indel=self.input()["indel"].path,
            input_snp=self.input()["snp"].path,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class new_Annovar1(Annovar1):
    def requires(self):
        return CombineVariants(sampleID=self.sampleID,
                               dry_run=self.dry_run)


class new_Annovar2(Annovar2):
    def requires(self):
        return new_Annovar1(sampleID=self.sampleID,
                            dry_run=self.dry_run)


class workflow(luigi.Task):
    x = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        samples_IDs = str(self.x).split(',')
        for i in samples_IDs:
            if NORMAL_SIG:
                if pfn(i, 'mt2_for') == NORMAL_SIG:
                    yield new_Annovar2(sampleID=i, dry_run=self.dry_run)
            else:
                yield new_Annovar2(sampleID=i, dry_run=self.dry_run)


if __name__ == '__main__':
    luigi.run()

    # python -m luigi --module SomaticPipelines_fast_version workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
