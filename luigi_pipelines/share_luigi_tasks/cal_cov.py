import luigi

class pre_analysis(luigi.Task):
    # todo: finish it add cal_cov part
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return GenerateSam_pair(infodict=self.infodict,
                                dry_run=self.dry_run)
    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.sam', '.bam'))

    def run(self):
        # fixme
        pass
