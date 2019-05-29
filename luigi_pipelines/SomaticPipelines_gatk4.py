from collections import defaultdict

import luigi

from luigi_pipelines.share_luigi_tasks import Annovar1, Annovar2
from luigi_pipelines.share_luigi_tasks.gatk4 import PrintReads
from . import *


#########somatic pipeline
class MuTect2_pair(luigi.Task):
    ####combine N & T
    sample_IDs = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        sampleIDs = str(self.sample_IDs).split(',')
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
        sampleIDs = str(self.sample_IDs).split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath,
                                                     PN=Project_ID,
                                                     PairN=pair_name) + '.mt2.bam'
        # todo :verify the suffix is right?
        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        prefix = self.output().path.rpartition('.bam')[0]
        input_normal = self.input()["normal"].path
        input_tumor = self.input()["tumor"].path
        if bed_file_path:
            extra_str = " --intervals %s" % bed_file_path
        else:
            extra_str = ''

        sampleIDs = str(self.sample_IDs).split(',')

        if pfn(sampleIDs[0], 'mt2_for') == NORMAL_SIG:
            normal_name = pfn(sampleIDs[0], 'sample_name')
            tumor_name = pfn(sampleIDs[1], 'sample_name')
        elif pfn(sampleIDs[0], 'mt2_for') == TUMOR_SIG:
            normal_name = pfn(sampleIDs[1], 'sample_name')
            tumor_name = pfn(sampleIDs[0], 'sample_name')
        else:
            raise Exception

        cmdline = "{gatk4} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_normal} -normal {N_name} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam {extra_str}".format(
            REF=REF_file_path,
            cosmic=cos_snp,
            db_snp=db_snp,
            input_tumor=input_tumor,
            input_normal=input_normal,
            gatk4=gatk_pro,
            N_name=normal_name,
            T_name=tumor_name,
            prefix=prefix,
            extra_str=extra_str)
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
        input1 = self.input()[0].path
        mt2_id = pfn(self.sample_NT, 'mt2_for')
        prefix = self.output().path.rpartition('.bam')[0]
        output_dir = self.output().path.rpartition('/')[0]

        if os.path.isdir(output_dir) != True:
            os.makedirs(output_dir)

        if bed_file_path:
            extra_str = " --intervals %s" % bed_file_path
        else:
            extra_str = ""

        if mt2_id == NORMAL_SIG:
            extra_str += ""
            # Normal only
        else:
            extra_str += " --tumor-lod-to-emit 4"
            # Tumor only

        cmdline = "{gatk4} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam".format(
            gatk4=gatk_pro,
            REF=REF_file_path,
            db_snp=db_snp,
            input_tumor=input1,
            prefix=prefix,
            T_name=pfn(self.sample_NT, 'sample_name'),
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
            import ipdb;
            ipdb.set_trace()
        for i in samples_IDs:
            yield new_Annovar2(sample_ID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
