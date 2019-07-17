import sys
from os.path import dirname, join

import luigi

sys.path.insert(0, dirname(dirname(__file__)))
sys.path.insert(0, dirname(__file__))
from parse_file_name import fileparser
from share_luigi_tasks import quality_assessment
from luigi_pipelines import run_cmd,valid_path
import click


project_root_path = dirname(dirname(__file__))
@click.group()
def cli():
    pass


@cli.command()
@click.argument('cmd', nargs=-1)
def run(cmd):
    luigi.run(cmdline_args=cmd)



@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_germline(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py workflow --tab {project_root_path}/test_set/germline/data_input.tsv --odir {odir} --analysis-type germline --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_germline_gatk4(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py workflow --tab {project_root_path}/test_set/germline/data_input.tsv --odir {odir} --analysis-type germline_gatk4 --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_somatic(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py workflow --tab {project_root_path}/test_set/somatic/data_input.tsv --odir {odir} --analysis-type somatic --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_somatic_gatk4(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py workflow --tab {project_root_path}/test_set/somatic/data_input.tsv --odir {odir} --analysis-type somatic_gatk4 --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)



class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="germline")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)

    def requires(self):
        df = fileparser(self.tab)
        antype = str(self.analysis_type).lower()
        tasks = {"germline_annovar2": [],
                 "somatic_annovar2": [],
                 "quality_accessment": [],
                 "somatic_pcgr": []}

        sample_info = []
        pair_info = []
        ############################################################
        if antype in ['germline', "germline_gatk4"]:
            if antype == "germline":
                from GermlinePipelines import new_Annovar2
            elif antype == "germline_gatk4":
                from GermlinePipelines_gatk4 import new_Annovar2
            else:
                raise Exception()

            for sample_name, sample_info in df.germline_pair().items():
                sample_info["odir"] = self.odir
                sample_info["log_path"] = self.log_path
                tasks["germline_annovar2"].append(new_Annovar2(infodict=sample_info,
                                                               dry_run=self.dry_run))

        elif antype in ['somatic', "somatic_gatk4"]:
            from post_pipelines_analysis.luigi_tasks.post_somatic import run_pcgr
            if antype == "somatic":
                from SomaticPipelines import new_Annovar2 as somatic_a2
                from GermlinePipelines import new_Annovar2 as germline_a2
            elif antype == "somatic_gatk4":
                from SomaticPipelines_gatk4 import new_Annovar2 as somatic_a2
                from GermlinePipelines_gatk4 import new_Annovar2 as germline_a2
            else:
                raise Exception("wrong somatic_gatk4")

            for combine_name, pair_info in df.get_full_info(self.odir, gettype='somatic').items():
                pair_info["odir"] = self.odir
                pair_info["source_name"] = pair_info["Normal"]["source_name"]
                pair_info["log_path"] = self.log_path
                pair_info["Normal"]["odir"] = self.odir
                pair_info["Normal"]["log_path"] = self.log_path
                pair_info["Tumor"]["odir"] = self.odir
                pair_info["Tumor"]["log_path"] = self.log_path

                tasks["somatic_annovar2"].append(somatic_a2(infodict=pair_info,
                                                            dry_run=self.dry_run))
                tasks["germline_annovar2"].append(germline_a2(infodict=pair_info["Normal"],
                                                              dry_run=self.dry_run))
                tasks["somatic_pcgr"].append(run_pcgr(infodict=pair_info,
                                                            dry_run=self.dry_run))
                # todo: adjust the post analysis module
        else:
            raise Exception

        if sample_info:
            # for germline (quality assessment)
            tasks["quality_accessment"].append(quality_assessment(tab_file=self.tab,
                                                                  odir=self.odir,
                                                                  infodict=sample_info,
                                                                  dry_run=self.dry_run,
                                                                  ))
            # quality assessment will generate cov.info and cov_summary.info
            # for two kinds of bam(sorted and recaled).
        if pair_info:
            # for somatic (quality assessment)
            tasks["quality_accessment"].append(quality_assessment(tab_file=self.tab,
                                                                  odir=self.odir,
                                                                  infodict=pair_info,
                                                                  dry_run=self.dry_run,
                                                                  ))

        return tasks

    def run(self):
        from post_analysis_tasks import germline_filter
        antype = str(self.analysis_type).lower()
        if 'germline' in antype:
            infile = self.input()["germline_annovar2"][0].path
            odir = join(dirname(infile), "germline_filter")
            germline_filter(dirname(infile),
                            odir,
                            self.tab)

        elif "somatic" in antype:
            pass
            # qa_files = self.input()["quality_accessment"]
            # csv_files = self.input()["somatic_annovar2"]
            # for combined_name, pair_info in df.get_full_info(self.odir, gettype='somatic').items():
            #     pair_info["log_path"] = self.log_path
            #     pair_info["odir"] = self.odir
            #     pair_info["source_name"] = pair_info["Normal"]["source_name"]
            #     pair_info["Normal"]["odir"] = self.odir
            #     pair_info["Normal"]["log_path"] = self.log_path
            #     pair_info["Tumor"]["odir"] = self.odir
            #     pair_info["Tumor"]["log_path"] = self.log_path
            #     luigi.build([preprocess_vcf(infodict=pair_info,
            #                           dry_run=self.dry_run)])


if __name__ == '__main__':
    cli()

# python3 luigi_pipelines/main.py workflow --tab test_set/germline/data_input.csv --odir test_set/germline_run_gatk4 --analysis-type germline --workers 5 --log-path test_set/germline_run_gatk4/cmd_log.txt

# python3 /home/liaoth/project/Whole_pipelines/luigi_pipelines/main.py workflow --tab /home/liaoth/project/Whole_pipelines/test_set/somatic/data_input.csv --odir /home/liaoth/project/Whole_pipelines/test_set/somatic_run --analysis-type somatic --workers 5 --log-path /home/liaoth/project/Whole_pipelines/test_set/somatic_run/cmd_log.txt
