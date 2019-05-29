import sys
from os.path import dirname

import luigi

sys.path.insert(0, dirname(dirname(__file__)))
from parse_file_name import fileparser


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    analysis_type = luigi.Parameter(default="germline")
    dry_run = luigi.BoolParameter()

    log_path = luigi.Parameter(default=None)

    df = fileparser(tab)

    def requires(self):
        if str(self.analysis_type).lower() == 'germline':
            from GermlinePipelines import new_Annovar2

            tasks = []
            for sample_name, sample_info in self.df.germline_pair().items():
                tasks.append(new_Annovar2(sample_dict=sample_info,
                                          dry_run=self.dry_run))
            return tasks

        elif str(self.analysis_type).lower() == 'somatic':
            from SomaticPipelines import new_Annovar2

            tasks = []
            for sample_name, sample_info in self.df.somatic_pair().items():
                tasks.append(new_Annovar2(sample_dict=sample_info,
                                          dry_run=self.dry_run))
            return tasks

        elif str(self.analysis_type).lower() == 'somatic_gemini':
            pass
