import sys
from os.path import dirname

import luigi

sys.path.insert(0, dirname(dirname(__file__)))
sys.path.insert(0, dirname(__file__))
from parse_file_name import fileparser
from share_luigi_tasks import quality_assessment


class main_entry(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="germline")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)

    def requires(self):
        df = fileparser(self.tab)
        antype = str(self.analysis_type).lower()
        tasks = []

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
                tasks.append(new_Annovar2(infodict=sample_info,
                                          dry_run=self.dry_run))


        elif antype in ['somatic', "somatic_gatk4"]:
            if antype == "somatic":
                from SomaticPipelines import new_Annovar2
            elif antype == "somatic_gatk4":
                from SomaticPipelines_gatk4 import new_Annovar2
            else:
                raise Exception("wrong somatic_gatk4")

            for combine_name, pair_info in df.get_full_info(self.odir, gettype='somatic').items():
                pair_info["log_path"] = self.log_path
                pair_info["Normal"]["odir"] = self.odir
                pair_info["Normal"]["log_path"] = self.log_path
                pair_info["Tumor"]["odir"] = self.odir
                pair_info["Tumor"]["log_path"] = self.log_path

                tasks.append(new_Annovar2(infodict=pair_info,
                                          dry_run=self.dry_run))

        # elif antype == "germline_pcgr":
        #     from GermlinePipelines_to_pcgr import gemini_part
        #     tasks = []
        #     for sample_name, sample_info in df.get_full_info(self.odir).items():
        #         sample_info["odir"] = self.odir
        #         sample_info["log_path"] = self.log_path
        #         tasks.append(gemini_part(infodict=sample_info,
        #                                  dry_run=self.dry_run))
        # elif antype == "somatic_pcgr":
        #     from SomaticPipelines_to_pcgr import gemini_part
        #
        #     for sample_name, pair_info in df.get_full_info(self.odir).items():
        #         sample_info["odir"] = self.odir
        #         sample_info["log_path"] = self.log_path
        #         tasks.append(gemini_part(infodict=sample_info,
        #                                  dry_run=self.dry_run))

        else:
            raise Exception

        if sample_info:
            # for germline
            tasks.append(quality_assessment(tab_file=self.tab,
                                            odir=self.odir,
                                            infodict=sample_info,
                                            dry_run=self.dry_run,
                                            ))
        if pair_info:
            # for somatic
            tasks.append(quality_assessment(tab_file=self.tab,
                                            odir=self.odir,
                                            infodict=pair_info,
                                            dry_run=self.dry_run,
                                            ))
        return tasks

    def run(self):
        pass


if __name__ == '__main__':
    luigi.run()

# python3 luigi_pipelines/main.py main_entry --tab test_set/germline/data_input.csv --odir test_set/germline_run_gatk4 --analysis-type germline --workers 5 --log-path test_set/germline_run_gatk4/cmd_log.txt

# python3 /home/liaoth/project/Whole_pipelines/luigi_pipelines/main.py main_entry --tab /home/liaoth/project/Whole_pipelines/test_set/somatic/data_input.csv --odir /home/liaoth/project/Whole_pipelines/test_set/somatic_run --analysis-type somatic --workers 5 --log-path /home/liaoth/project/Whole_pipelines/test_set/somatic_run/cmd_log.txt

