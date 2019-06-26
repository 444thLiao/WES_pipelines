import sys
from os.path import dirname, join

sys.path.insert(0, dirname(dirname(dirname(__file__))))
import luigi
from luigi_pipelines import run_cmd, valid_path, config
from luigi_pipelines.share_luigi_tasks import BaseRecalibrator, base_luigi_task
from luigi_pipelines.SomaticPipelines_gatk4 import new_Annovar2 as somatic_a2
from luigi_pipelines.SomaticPipelines_gatk4 import new_Annovar1 as somatic_a1
from GermlinePipelines_gatk4 import new_Annovar2 as germline_a2

project_root_path = dirname(dirname(dirname(__file__)))


class add_cov(base_luigi_task):
    # infodict is a pair input. not one of them
    def requires(self):
        tasks = {}
        tasks["normal"] = BaseRecalibrator(infodict=self.infodict["Normal"],
                                           dry_run=self.dry_run)
        tasks["tumor"] = BaseRecalibrator(infodict=self.infodict["Tumor"],
                                          dry_run=self.dry_run)
        tasks["annovar"] = somatic_a2(infodict=self.infodict,
                                      dry_run=self.dry_run)
        # this show include pair and tumor single output csv.
        return tasks

    def output(self):
        ofiles = []
        for i in self.input()["annovar"]:
            ofiles.append(i.path.replace('.csv',
                                         '_added_cov.csv'))
        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        normal_bam = self.input()["normal"].path
        tumor_bam = self.input()["tumor"].path
        input_csv_files = self.input()["annovar"]

        for input_csv in input_csv_files:
            input_csv = input_csv.path
            output_csv = input_csv.replace('.csv',
                                           '_added_cov.csv')
            cmdline = f"python3 {project_root_path}/api/add_per_info_into_csv.py -i {input_csv} -o {output_csv} -tb {tumor_bam} -nb {normal_bam}"
            run_cmd(cmdline, log_file=self.get_log_path(), dry_run=self.dry_run)


class perform_filters(base_luigi_task):

    def requires(self):
        tasks = {}
        tasks["pair&singleT"] = add_cov(infodict=self.infodict,
                                        dry_run=self.dry_run)
        tasks["normal_germline"] = germline_a2(infodict=self.infodict["Normal"],
                                               dry_run=self.dry_run)
        return tasks

    def output(self):
        odir = self.infodict["odir"]
        source_name = self.infodict["source_name"]
        filter_path = config.filtered_dir.format(path=odir,
                                                 PairN=source_name)
        valid_path(filter_path, check_odir=1)
        opath1 = join(filter_path, source_name + '_except_AF_depth.csv')
        opath2 = join(filter_path, source_name + '_except_AF_depth_PASS.csv')
        return [luigi.LocalTarget(opath1), luigi.LocalTarget(opath2)]

    def run(self):
        from post_pipelines_analysis.filter_pipelines import filter_pipelines2

        normal_germline_csv = self.input()["normal_germline"][0].path
        pair_somatic_csv = self.input()["pair&singleT"][0].path
        tumor_somatic_csv = self.input()["pair&singleT"][1].path
        # the order of its is defined at `share_luigi_tasks.annovar_tasks`
        if not self.dry_run:
            filter_pipelines2(normal_germline=normal_germline_csv,
                              tumor_somatic=tumor_somatic_csv,
                              pair_somatic=pair_somatic_csv,
                              pp=[3, 4, 5, 6],
                              output_path=self.output()[0].path)
            filter_pipelines2(normal_germline=normal_germline_csv,
                              tumor_somatic=tumor_somatic_csv,
                              pair_somatic=pair_somatic_csv,
                              pp=[3, 4, 6],
                              output_path=self.output()[1].path)


class preprocess_vcf(base_luigi_task):

    def requires(self):
        tasks = {}
        tasks["filtered"] = perform_filters(infodict=self.infodict,
                                            dry_run=self.dry_run)
        tasks["pair_vcf"] = somatic_a1(infodict=self.infodict,
                                       dry_run=self.dry_run,
                                       mode='pair')
        tasks["tumor_vcf"] = somatic_a1(infodict=self.infodict["Tumor"],
                                        dry_run=self.dry_run,
                                        mode='single')
        return tasks

    def output(self):
        odir = self.infodict["odir"]
        source_name = self.infodict["source_name"]
        filter_path = config.filtered_dir.format(path=odir,
                                                 PairN=source_name)
        valid_path(filter_path, check_odir=1)

        final_variants = join(filter_path, source_name + '_final.vcf')
        return luigi.LocalTarget(final_variants)

    def run(self):
        from post_pipelines_analysis.extracted_pos_from_vcf import merge_two_vcf
        from special_fun.csv2bed import csv2bed
        pair_vcf = self.input()["pair_vcf"].path.replace('.av','.vcf')
        tumor_single_vcf = self.input()["tumor_vcf"].path.replace('.av','.vcf')
        filtered_csv = self.input()["filtered"][1].path

        if not self.dry_run:
            filtered_bed = csv2bed(filtered_csv,
                                   filtered_csv.replace('.csv', '.bed'))
            merge_two_vcf(pair_vcf,
                          tumor_single_vcf,
                          filtered_bed,
                          self.output().path,
                          log_file=self.get_log_path())

        cmdline = "{vt} decompose -s {input_vcf} | {vt} normalize -r {REF_path} - > {output_vcf}".format(vt=config.vt_pro,
                                                                                                         input_vcf=self.output().path,
                                                                                                         REF_path=config.REF_file_path,
                                                                                                         output_vcf=self.output().path.replace('.vcf',
                                                                                                                                               '.vt.vcf'))
        run_cmd(cmdline, dry_run=self.dry_run, log_file=self.get_log_path())



class run_pcgr(base_luigi_task):
    def requires(self):
        return preprocess_vcf(infodict=self.infodict,
                              dry_run=self.dry_run)

    def output(self):
        indir = dirname(self.input().path)
        odir = join(indir,"pcgr_output")
        valid_path(odir,check_odir=1)
        source_name = self.infodict["source_name"]
        ofile = join(odir,"{}.pcgr_acmg.grch37.html".format(source_name))
        # todo: convert ref_path to grch37??
        return luigi.LocalTarget(ofile)

    def run(self):
        input_vcf = self.input().path.replace('.vcf', '.vt.vcf')
        output_dir = dirname(self.output().path)
        source_name = self.infodict["source_name"]
        valid_path(output_dir, check_odir=1)
        # todo: convert ref_path to grch37??
        cmdline = "source activate pcgr; python3 {pcgr_dir}/pcgr.py --input_vcf {input_vcf} {pcgr_dir} {output_dir} grch37 {toml_config} {source_name} --no-docker --force_overwrite".format(
            pcgr_dir=config.pcgr_dir,
            input_vcf=input_vcf,
            output_dir=output_dir,
            toml_config=config.pcgr_toml_file,
            source_name=source_name)
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())

class packitup(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):

        pass

    def output(self):
        pass

    # def run(self):
    #     somatic_pp_output_root_dir = str(self.odir)
    #     summary_p = join(somatic_pp_output_root_dir,
    #                      "summary_results")
    #     somatic_csv = join(summary_p, 'somatic/')
    #     germline_csv = join(summary_p, 'germline/')
    #     vcf_path = join(summary_p, 'vcf_storge/')
    #     # info_summary_path = join(summary_p, 'summary_info/')
    #     # # relative output path
    #     # final_csv = join(summary_p, 'final_output/')
    #     # filtered_csv = join(summary_p, 'filtered_file/')
    #     # pcgr_output = join(summary_p, 'pcgr/')
    #     # summary_depth = join(summary_p, 'info_summary/')
    #     # pack_result = join(summary_p, 'packet_result')
    #
    #     cmdline = f'cp {somatic_pp_output_root_dir}/*_result/*_somatic/*anno.csv {somatic_csv}; \
    #     cp {somatic_pp_output_root_dir}/*_result/N*/*anno.csv {germline_csv}; \
    #      cp {somatic_pp_output_root_dir}/*_result/*_somatic/*.mt2.vcf {vcf_path}; \
    #        '
    #     run_cmd(cmdline, dry_run=self.dry_run)
