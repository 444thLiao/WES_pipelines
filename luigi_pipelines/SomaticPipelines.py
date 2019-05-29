from collections import defaultdict

import luigi

from luigi_pipelines.share_luigi_tasks import PrintReads, run_cmd, valid_path, Annovar1, Annovar2
from main import *


#########somatic pipeline
class MuTect2_pair(luigi.Task):
    ####combine N & T
    sample_IDs = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        sampleIDs = self.sample_IDs.split(',')
        if pfn(sampleIDs[0], 'mt2_for') == NORMAL_SIG:
            normal_id = sampleIDs[0]
            tumor_id = sampleIDs[1]
        elif pfn(sampleIDs[0], 'mt2_for') == TUMOR_SIG:
            normal_id = sampleIDs[1]
            tumor_id = sampleIDs[0]
        else:
            raise Exception

        required_tasks = {}
        required_tasks["normal"] = PrintReads(sampleID=normal_id,
                                              dry_run=self.dry_run)
        required_tasks["tumor"] = PrintReads(sampleID=tumor_id,
                                             dry_run=self.dry_run)
        return required_tasks

    def output(self):
        # todo: simplify
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath,
                                                     PN=Project_ID,
                                                     PairN=pair_name) + '.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_normal = self.input()["normal"].path
        input_tumor = self.input()["tumor"].path

        prefix = self.output().path.rpartition('.bam')[0]
        cmdline = '''java -Xmx10g -jar {gatk} -T MuTect2 --allSitePLs -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:normal {input_normal} --input_file:tumor {input_tumor} --out {prefix}.vcf --bamOutput {output_f --log_to_file {prefix}.log'''.format(
            gatk=gatkv36_path,
            REF=REF_file_path,
            cosmic=cos_snp,
            db_snp=db_snp,
            input_tumor=input_tumor,
            input_normal=input_normal,
            prefix=prefix,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class MuTect2_single(luigi.Task):
    sample_NT = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return PrintReads(sampleID=self.sample_NT,
                          dry_run=self.dry_run)

    def output(self):

        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')

        output_path = somatic_single_output_fmt.format(path=base_outpath,
                                                       PN=Project_ID,
                                                       SN=sample_name) + '.mt2.bam'

        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_f = self.input().path
        mt2_id = pfn(self.sample_NT, 'mt2_for')
        prefix = self.output().path.rpartition('.bam')[0]

        if mt2_id == NORMAL_SIG:
            # Normal only
            extra_str = ''
        else:
            # Tumor only
            extra_str = ' --tumor_lod 4'
        # both normal and tumor sample use input_file:tumor as parameter
        cmdline = '''java -Xmx10g -jar {gatk} -T MuTect2 --allSitePLs --artifact_detection_mode -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_f} --out {prefix}.vcf --bamOutput {output_f} --log_to_file {prefix}.log {extra_str}'''.format(
            gatk=gatkv36_path,
            REF=REF_file_path,
            cosmic=cos_snp,
            db_snp=db_snp,
            input_f=input_f,
            output_f=self.output().path,
            prefix=prefix,
            extra_str=extra_str)
        run_cmd(cmdline, dry_run=self.dry_run)


class new_Annovar1(Annovar1):

    def requires(self):
        if self.sample_ID in pair_bucket:
            pair_value = pair_bucket[self.sample_ID]
            return [MuTect2_pair(sample_IDs=','.join(pair_value))]
        else:
            return [MuTect2_single(sample_NT=self.sample_ID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.bam')[0] + '.merged.av')


class new_Annovar2(Annovar2):
    def requires(self):
        return new_Annovar1(sampleID=self.sampleID,
                            dry_run=self.dry_run)


class workflow(luigi.Task):
    x = luigi.Parameter()

    def requires(self):
        samples_IDs = str(self.x).split(',')

        pair_bucket = defaultdict(list)
        for _x in samples_IDs:
            pair_bucket[pfn(_x,
                            'pair_name')].append(_x)
        global pair_bucket
        ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}

        samples_IDs += [_x for _x in pair_bucket.keys()]
        for i in samples_IDs:
            yield new_Annovar2(sample_ID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
