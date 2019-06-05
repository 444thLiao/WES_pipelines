import luigi

from luigi_pipelines import config
from toolkit import run_cmd, valid_path
from .basic_tasks import base_luigi_task

#########14
class Annovar1(base_luigi_task):

    def requires(self):
        "abstract for implement"
        pass

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vcf','.av'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        cmdline = "%s/convert2annovar.pl %s --includeinfo -format vcf4 > %s" % (
            config.annovar_pro,
            self.input().path,
            self.output().path)
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

class Annovar2(base_luigi_task):

    def requires(self):
        "abstract for implement"
        # must be a list (even contain one element)
        pass
        # return Annovar1(sample_dict=self.sample_dict,
        #                 dry_run=self.dry_run)

    def output(self):
        if type(self.input()) == dict:
            ofiles = [luigi.LocalTarget(
                _input.path.replace('.av',
                                    '.anno.%s_multianno.csv' % config.genome_version))
                for _input in self.input().values()]
        elif type(self.input()) == list:
            ofiles = [luigi.LocalTarget(
                _input.path.replace('.av',
                                    '.anno.%s_multianno.csv' % config.genome_version))
                for _input in self.input()]
        else:
            return [luigi.LocalTarget(
                self.input().path.replace('.av',
                                          '.anno.%s_multianno.csv' % config.genome_version))]
        return ofiles

    def run(self):
        if type(self.input()) == dict:
            input_list = list(self.input().values())
        elif type(self.input()) == list:
            input_list = list(self.input())
        else:
            input_list = [self.input()]

        for _output, _input in zip(self.output(), input_list):
            valid_path(_output.path, check_ofile=1)
            prefix = _input.path.replace('.av','')
            cmdline = "{annovar_dir}/table_annovar.pl {input_f} {annovar_db} -buildver {genome_version} -protocol {db_names} -operation g,r,r,f,f,f,f,f,f -nastring . --remove --otherinfo --csvout --thread {annovar_thread} --outfile {output_f} --argument '-exonicsplicing -splicing 25',,,,,,,,".format(
                annovar_dir=config.annovar_pro,
                input_f=_input.path,
                annovar_db=config.annovar_db,
                genome_version=config.genome_version,
                db_names=config.db_names,
                annovar_thread=config.annovar_thread,
                output_f=prefix + '.anno')
            run_cmd(cmdline,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path())
            if self.dry_run:
                run_cmd("touch %s" % _output.path, dry_run=False)