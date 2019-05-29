import luigi

from luigi_pipelines import run_cmd, valid_path, config
from luigi_pipelines.share_luigi_tasks import PrintReads, Annovar1, Annovar2


#########somatic pipeline
class MuTect2_pair(luigi.Task):
    ####combine N & T
    infodict_N = luigi.DictParameter()
    infodict_T = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        required_tasks = {}
        required_tasks["normal"] = PrintReads(infodict=self.infodict_N,
                                              dry_run=self.dry_run)
        required_tasks["tumor"] = PrintReads(infodict=self.infodict_T,
                                             dry_run=self.dry_run)
        return required_tasks

    def output(self):
        project_name = self.infodict.get("project_name", "")
        odir = self.infodict.get("odir", "")

        source_name = self.infodict.get("source_name", "")
        output_path = config.somatic_pair_output_fmt.format(
            path=odir,
            PN=project_name,
            PairN=source_name) + '.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_normal = self.input()["normal"].path
        input_tumor = self.input()["tumor"].path

        prefix = self.output().path.rpartition('.bam')[0]

        cmdline = '''java -Xmx10g -jar {gatk} -T MuTect2 --allSitePLs -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:normal {input_normal} --input_file:tumor {input_tumor} --out {prefix}.vcf --bamOutput {output_f --log_to_file {prefix}.log'''.format(
            gatk=config.gatkv36_path,
            REF=config.REF_file_path,
            cosmic=config.cos_snp,
            db_snp=config.db_snp,
            input_tumor=input_tumor,
            input_normal=input_normal,
            prefix=prefix,
            output_f=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class MuTect2_single(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return PrintReads(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        sample_name = self.infodict.get("SampleID", '')
        project_name = self.infodict.get("project_name", "")
        odir = self.infodict.get("odir", "")

        output_path = config.somatic_single_output_fmt.format(
            path=odir,
            PN=project_name,
            SN=sample_name) + '.mt2.bam'

        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_f = self.input().path

        somatic_type = self.infodict.get("Somatic")
        if somatic_type == "N":
            # Normal only
            extra_str = ''
        elif somatic_type == "T":
            # Tumor only
            extra_str = ' --tumor_lod 4'
        else:
            raise Exception("Unknown values of Somatic columns (like '%s' )" % somatic_type)
        # both normal and tumor sample use input_file:tumor as parameter
        cmdline = '''java -Xmx10g -jar {gatk} -T MuTect2 --allSitePLs --artifact_detection_mode -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_f} --out {prefix}.vcf --bamOutput {output_f} --log_to_file {prefix}.log {extra_str}'''.format(
            gatk=config.gatkv36_path,
            REF=config.REF_file_path,
            cosmic=config.cos_snp,
            db_snp=config.db_snp,
            input_f=input_f,
            output_f=self.output().path,
            prefix=input_f.replace('.bam', ''),
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
#             pair_bucket[pfn(_x,
#                             'pair_name')].append(_x)
#         global pair_bucket
#         ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}
#
#         samples_IDs += [_x for _x in pair_bucket.keys()]
#         for i in samples_IDs:
#             yield new_Annovar2(sample_ID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
