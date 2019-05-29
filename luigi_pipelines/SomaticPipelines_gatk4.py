import luigi

from luigi_pipelines import config, run_cmd, valid_path
from luigi_pipelines.SomaticPipelines import MuTect2_pair, MuTect2_single
from luigi_pipelines.share_luigi_tasks import Annovar1, Annovar2
from luigi_pipelines.share_luigi_tasks.gatk4 import PrintReads


#########somatic pipeline
class MuTect2_pair(MuTect2_pair):
    ####combine N & T

    def requires(self):
        required_tasks = {}
        required_tasks["normal"] = PrintReads(infodict=self.infodict_N,
                                              dry_run=self.dry_run)
        required_tasks["tumor"] = PrintReads(infodict=self.infodict_T,
                                             dry_run=self.dry_run)
        return required_tasks

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        prefix = self.output().path.replace('.bam', "")
        input_normal = self.input()["normal"].path
        input_tumor = self.input()["tumor"].path
        if config.bed_file_path:
            extra_str = " --intervals %s" % config.bed_file_path
        else:
            extra_str = ''

        normal_name = self.infodict_N["SampleID"]
        tumor_name = self.infodict_T["SampleID"]

        cmdline = "{gatk4} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_normal} -normal {N_name} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam {extra_str}".format(
            REF=config.REF_file_path,
            cosmic=config.cos_snp,
            db_snp=config.db_snp,
            input_tumor=input_tumor,
            input_normal=input_normal,
            gatk4=config.gatk_pro,
            N_name=normal_name,
            T_name=tumor_name,
            prefix=prefix,
            extra_str=extra_str)
        run_cmd(cmdline, dry_run=self.dry_run)


class MuTect2_single(MuTect2_single):
    def requires(self):
        return PrintReads(infodict=self.infodict,
                          dry_run=self.dry_run)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_f = self.input().path
        sample_name = self.infodict.get("SampleID", '')
        if config.bed_file_path:
            extra_str = " --intervals %s" % config.bed_file_path
        else:
            extra_str = ""

        somatic_type = self.infodict.get("Somatic")
        if somatic_type == "N":
            extra_str += ""
            # Normal only
            extra_str = ''
        elif somatic_type == "T":
            # Tumor only
            extra_str = ' --tumor_lod 4'
        else:
            raise Exception("Unknown values of Somatic columns (like '%s' )" % somatic_type)

        cmdline = "{gatk4} Mutect2 --java-options '-Xmx20g' --native-pair-hmm-threads 20 --reference {REF} -I {input_tumor} -tumor {T_name} --dbsnp {db_snp} --seconds-between-progress-updates 60 --all-site-pls -stand-call-conf 10 -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A BaseQuality -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts --all-site-pls true --output {prefix}.vcf -bamout {prefix}.bam".format(
            gatk4=config.gatk_pro,
            REF=config.REF_file_path,
            db_snp=config.db_snp,
            input_tumor=input_f,
            prefix=input_f.replace('.bam', ''),
            T_name=sample_name,
            extra_str=extra_str)
        run_cmd(cmdline, dry_run=self.dry_run)


class new_Annovar1(Annovar1):

    def requires(self):
        normal_dict = self.infodict["Normal"]
        tumor_dict = self.infodict["Tumor"]

        tasks = {}
        tasks["pair"] = MuTect2_pair(infodict_N=normal_dict,
                                     infodict_T=tumor_dict,
                                     dry_run=self.dry_run)
        tasks["single_N"] = MuTect2_single(infodict=normal_dict,
                                           dry_run=self.dry_run)
        tasks["single_T"] = MuTect2_single(infodict=tumor_dict,
                                           dry_run=self.dry_run)
        return tasks

    def output(self):
        tasks = []
        for localtarget in self.input():
            tasks.append(luigi.LocalTarget(localtarget.path.rpartition('.bam')[0] + '.merged.av'))
        return tasks


class new_Annovar2(Annovar2):
    def requires(self):
        return new_Annovar1(infodict=self.infodict,
                            dry_run=self.dry_run)

#
# class workflow(luigi.Task):
#     x = luigi.Parameter()
#
#     def requires(self):
#         samples_IDs = str(self.x).split(',')
#
#         pair_bucket = defaultdict(list)
#         for _x in samples_IDs:
#             pair_bucket[pfn(_x, 'pair_name')].append(_x)
#         adjust_multiple = []
#         for each in pair_bucket.keys():
#             if len(pair_bucket[each]) > 2:
#                 tmp = pair_bucket[each]
#                 only_normal = [_ for _ in tmp if pfn(_, 'mt2_for') == NORMAL_SIG][0]
#                 for _each in tmp:
#                     if pfn(_each, 'mt2_for') == TUMOR_SIG and pfn(_each, 'sample_name').replace(TUMOR_SIG, '') != each:
#                         adjust_multiple.append((pfn(_each, 'sample_name').replace(TUMOR_SIG, ''), [only_normal, _each]))
#                     elif pfn(_each, 'mt2_for') == TUMOR_SIG and pfn(_each, 'sample_name').replace(TUMOR_SIG,
#                                                                                                   '') == each:
#                         adjust_multiple.append((each, [only_normal, _each]))
#         pair_bucket.update(dict(adjust_multiple))
#         global pair_bucket
#         ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}
#
#         samples_IDs += [_x for _x in pair_bucket.keys()]
#
#         if debug_:
#             import ipdb;
#             ipdb.set_trace()
#         for i in samples_IDs:
#             yield new_Annovar2(sample_ID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
