import sys
from os.path import dirname

import luigi

sys.path.insert(0, dirname(dirname(__file__)))
from parse_file_name import fileparser


class main_entry(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="germline")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)

    def requires(self):
        df = fileparser(self.tab)

        if str(self.analysis_type).lower() == 'germline':
            from GermlinePipelines import new_Annovar2

            tasks = []
            for sample_name, sample_info in df.germline_pair().items():
                sample_info["odir"] = self.odir
                sample_info["log_path"] = self.log_path
                tasks.append(new_Annovar2(sample_dict=sample_info,
                                          dry_run=self.dry_run))
            return tasks

        elif str(self.analysis_type).lower() == 'somatic':
            from SomaticPipelines import new_Annovar2

            tasks = []
            for sample_name, sample_info in df.somatic_pair().items():
                tasks.append(new_Annovar2(sample_dict=sample_info,
                                          dry_run=self.dry_run))
            return tasks

        elif str(self.analysis_type).lower() == 'somatic_gemini':
            pass
if __name__ == '__main__':
    luigi.run()

# python3 main_entry --tab /home/liaoth/project/ZHJ_WES/data_input.csv --odir /home/liaoth/project/ZHJ_WES/output --analysis-type germline --dry-run