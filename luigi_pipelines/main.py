import sys
from os.path import dirname

import luigi

sys.path.insert(0, dirname(dirname(__file__)))
sys.path.insert(0, dirname(__file__))
from parse_file_name import fileparser


class main_entry(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="germline")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)

    def requires(self):
        df = fileparser(self.tab)
        antype = str(self.analysis_type).lower()
        if antype in ['germline', "germline_gatk4"]:
            if antype == "germline":
                from GermlinePipelines import new_Annovar2
            elif antype == "germline_gatk4":
                from GermlinePipelines_gatk4 import new_Annovar2
            else:
                raise Exception()
            tasks = []
            for sample_name, sample_info in df.get_output_file_path(self.odir).items():
                sample_info["odir"] = self.odir
                sample_info["log_path"] = self.log_path
                tasks.append(new_Annovar2(infodict=sample_info,
                                          dry_run=self.dry_run))
            return tasks

        elif antype in ['somatic', "somatic_gatk4"]:
            if antype == "germline":
                from SomaticPipelines import new_Annovar2
            elif antype == "germline_gatk4":
                from SomaticPipelines_gatk4 import new_Annovar2
            else:
                raise Exception()

            tasks = []
            for combine_name, pair_info in df.somatic_pair().items():
                pair_info["odir"] = self.odir
                pair_info["log_path"] = self.log_path
                pair_info["Normal"]["odir"] = self.odir
                pair_info["Normal"]["log_path"] = self.log_path
                pair_info["Tumor"]["odir"] = self.odir
                pair_info["Tumor"]["log_path"] = self.log_path
                tasks.append(new_Annovar2(infodict=pair_info,
                                          dry_run=self.dry_run))
            return tasks
        elif antype == "germline_gemini":
            from GermlinePipelines_to_gemini import gemini_part
            tasks = []
            for sample_name, sample_info in df.get_output_file_path(self.odir).items():
                sample_info["odir"] = self.odir
                sample_info["log_path"] = self.log_path
                tasks.append(gemini_part(infodict=sample_info,
                                             dry_run=self.dry_run))
            print(tasks)
            return tasks
        elif antype == "somatic_gemini":
            from SomaticPipelines_to_gemini import gemini_part
            tasks = []
            for sample_name, sample_info in df.get_output_file_path(self.odir).items():
                sample_info["odir"] = self.odir
                sample_info["log_path"] = self.log_path
                tasks.append(gemini_part(infodict=sample_info,
                                             dry_run=self.dry_run))

            return tasks
        else:
            raise Exception



if __name__ == '__main__':
    luigi.run()

# python3 main_entry --tab /home/liaoth/project/ZHJ_WES/data_input.csv --odir /home/liaoth/project/ZHJ_WES/output --analysis-type germline --dry-run
