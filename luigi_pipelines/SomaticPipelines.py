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
        project_name = self.infodict_N["project_name"]
        odir = self.infodict_N["odir"]

        source_name = self.infodict_N["source_name"]
        output_path = config.somatic_pair_output_fmt.format(
            path=odir,
            PN=project_name,
            PairN=source_name) + '.mt2.vcf'
        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        prefix = self.output().path.replace('.vcf', "")
        input_normal = self.input()["normal"].path
        input_tumor = self.input()["tumor"].path


        cmdline = '''java {java_option} -jar {gatk} -T MuTect2 --allSitePLs -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:normal {input_normal} --input_file:tumor {input_tumor} --out {output_f} --bamOutput {prefix}.bam --log_to_file {prefix}.log'''.format(
            gatk=config.gatkv36_path,
            java_option=config.java_option,
            REF=config.REF_file_path,
            cosmic=config.cos_snp,
            db_snp=config.db_snp,
            input_tumor=input_tumor,
            input_normal=input_normal,
            prefix=prefix,
            output_f=self.output().path)
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.infodict_N.get("log_path",None))
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

class MuTect2_single(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return PrintReads(infodict=self.infodict,
                          dry_run=self.dry_run)

    def output(self):
        sample_name = self.infodict["SampleID"]
        project_name = self.infodict["project_name"]
        odir = self.infodict["odir"]

        output_path = config.somatic_single_output_fmt.format(
            path=odir,
            PN=project_name,
            SN=sample_name) + '.mt2.vcf'

        return luigi.LocalTarget(output_path)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        input_f = self.input().path

        somatic_type = self.infodict["Somatic"]
        if somatic_type == "N":
            # Normal only
            extra_str = ''
        elif somatic_type == "T":
            # Tumor only
            extra_str = ' --tumor_lod 4'
        else:
            raise Exception("Unknown values of Somatic columns (like '%s' )" % somatic_type)
        # both normal and tumor sample use input_file:tumor as parameter
        cmdline = '''java {java_option} -jar {gatk} -T MuTect2 --allSitePLs --artifact_detection_mode -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_f} --out {output_f} --bamOutput {prefix}.bam --log_to_file {prefix}.log {extra_str}'''.format(
            gatk=config.gatkv36_path,
            java_option=config.java_option,
            REF=config.REF_file_path,
            cosmic=config.cos_snp,
            db_snp=config.db_snp,
            input_f=input_f,
            output_f=self.output().path,
            prefix=self.output().path.replace('.vcf', ''),
            extra_str=extra_str)
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.infodict.get("log_path",None))
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

class new_Annovar1(Annovar1):
    mode = luigi.Parameter()

    def requires(self):
        if self.mode == 'pair':

            normal_dict = self.infodict["Normal"]
            tumor_dict = self.infodict["Tumor"]

            return MuTect2_pair(infodict_N=normal_dict,
                                infodict_T=tumor_dict,
                                dry_run=self.dry_run)
        elif self.mode == 'single':

            return MuTect2_single(infodict=self.infodict,
                                  dry_run=self.dry_run)
        else:
            raise Exception("wrong input")

class new_Annovar2(Annovar2):

    def requires(self):
        tasks = {}
        tasks["pair"] = new_Annovar1(infodict=self.infodict,
                                     dry_run=self.dry_run,
                                     mode='pair')
        tasks["single_N"] = new_Annovar1(infodict=self.infodict["Normal"],
                                         dry_run=self.dry_run,
                                         mode='single')
        tasks["single_T"] = new_Annovar1(infodict=self.infodict["Tumor"],
                                         dry_run=self.dry_run,
                                         mode='single')
        return tasks


# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
