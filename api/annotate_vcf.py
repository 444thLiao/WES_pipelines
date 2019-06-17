"""
It also an api demo for someone who want to use part of these module.
"""
import os

import sys
from os.path import dirname

import luigi
sys.path.insert(0, dirname(dirname(__file__)))
from api import Annovar2, Annovar1


class input_file(luigi.ExternalTask):
    """
    given a input vcf. suffix must be vcf
    """
    vcf_file = luigi.Parameter()

    def output(self):
        vcf_path = os.path.abspath(self.vcf_file)
        if not vcf_path.endswith(".vcf"):
            raise Exception("input vcf file must end with .vcf")
        return luigi.LocalTarget(vcf_path)


class Annovar1(Annovar1):
    vcf_file = luigi.Parameter(default=None)
    output_dir = luigi.Parameter(default=None)
    def requires(self):
        return input_file(vcf_file=self.vcf_file)

    def output(self):
        if self.output_dir is None:
            ofile = self.input().path.replace('.vcf','.av')
        else:
            ofile = os.path.join(os.path.abspath(self.output_dir),
                                 os.path.basename(self.input().path).replace('.vcf', '.av'))
        return luigi.LocalTarget(ofile)

class Annovar2(Annovar2):
    log_file = luigi.Parameter(default=None)
    vcf_file = luigi.Parameter(default=None)
    output_dir = luigi.Parameter(default=None)
    def requires(self):
        if self.log_file is not None:
            log_path = os.path.abspath(self.log_file)
            info_dict = {"log_path": log_path}
        else:
            info_dict = {}
        return Annovar1(infodict=info_dict,
                        vcf_file=self.vcf_file,
                        dry_run=self.dry_run)


if __name__ == '__main__':
    luigi.run()

    # python3 /home/liaoth/data2/project/Whole_pipelines/api/annotate_vcf.py Annovar2 --vcf-file /home/liaoth/Desktop/someone.mt2.vcf --output-dir /home/liaoth/Desktop/
