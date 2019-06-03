import luigi

from luigi_pipelines import config
from toolkit import run_cmd, valid_path


#########14
class Annovar1(luigi.Task):
    infodict = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        "abstract for implement"
        pass

    def output(self):
        return luigi.LocalTarget(self.input().path.rpartition('.bam')[0] + '.merged.av')

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        cmdline = "%s/convert2annovar.pl %s --includeinfo -format vcf4 > %s" % (
            config.annovar_pro,
            self.input().path,
            self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class Annovar2(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        "abstract for implement"
        # must be a list (even contain one element)
        pass
        # return Annovar1(sample_dict=self.sample_dict,
        #                 dry_run=self.dry_run)

    def output(self):
        if type(self.input()) == dict:
            ofiles = [luigi.LocalTarget(
                _input.path.replace('.merged.av',
                                    '.merged.anno.%s_multianno.csv' % config.genome_version))
                for _input in self.input().values()]
        elif type(self.input()) == list:
            ofiles = [luigi.LocalTarget(
                _input.path.replace('.merged.av',
                                    '.merged.anno.%s_multianno.csv' % config.genome_version))
                for _input in self.input()]
        else:
            return [luigi.LocalTarget(
                self.input().path.replace('.merged.av',
                                          '.merged.anno.%s_multianno.csv' % config.genome_version))]
        return ofiles

    def run(self):
        for _output, _input in zip(self.output(), self.input()):
            valid_path(_output.path, check_ofile=1)
            prefix = _input.path.rpartition('.merged.av')[0]
            cmdline = "{annovar_dir}/table_annovar.pl {input_f} {annovar_db} -buildver {genome_version} -protocol {db_names} -operation g,r,r,f,f,f,f,f,f -nastring . --remove --otherinfo --csvout --thread {annovar_thread} --outfile {output_f} --argument '-exonicsplicing -splicing 25',,,,,,,,".format(
                annovar_dir=config.annovar_pro,
                input_f=_input.path,
                annovar_db=config.annovar_db,
                genome_version=config.genome_version,
                db_names=config.db_names,
                annovar_thread=config.annovar_thread,
                output_f=prefix + '.merged.anno')
            run_cmd(cmdline, dry_run=self.dry_run)
