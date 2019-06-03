import luigi

from pre_pipelines_analysis.cal_Cov_script_version import bed2info
from special_fun import Add_cov_ino_in_vcf as P_vcf
from special_fun.vcf_2_bed import vcf2bed
from luigi_pipelines import config, run_cmd


class luigi_vcf2bed(luigi.Task):
    "convert vcf to bed, summarized coverage info from this bed"
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        pass
        # return CombineVariants(infodict=self.infodict, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vcf',
                                                           '.bed'))

    def run(self):
        vcf2bed(self.input().path,
                self.output().path)


class luigi_bed2info(luigi.Task):
    "convert vcf to bed, summarized coverage info from this bed"
    infodict = luigi.DictParameter()
    ref_file = luigi.Parameter(default='')
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        pass
        # return [PrintReads(infodict=self.infodict, dry_run=self.dry_run),
        #         convert_vcf2bed(infodict=self.infodict, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace('.bed',
                                                              '.info'))

    def run(self):
        bed2info(self.input()[0].path,
                 self.output().path,
                 self.input()[1].path,
                 config.REF_file_path if not self.ref_file else self.ref_file)


#########14
class Add_cov_infos(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        "abstract for implement"
        pass
        # return [CombineVariants(infodict=self.infodict,dry_run=self.dry_run),
        #         luigi_bed2info(infodict=self.infodict,dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.vcf',
                                                               '.with_cov.vcf'))

    def run(self):
        merged_vcf = self.input()[0].path
        cov_info = self.input()[1].path

        P_vcf.Add_in_vcf_SO(infofile=cov_info,
                            vcf_path=merged_vcf,
                            output_vcf=self.output().path,)


#########15
class vt_part(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return Add_cov_infos(infodict=self.infodict, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.with_cov.vcf',
                                                           '.vt.vcf'))

    def run(self):
        cmdline = "{vt} decompose -s {input_vcf} | {vt} normalize -r {REF} - > {vt_vcf}".format(
            vt=config.vt_pro,
            input_vcf=self.input().path,
            REF=config.REF_file_path,
            vt_vcf=self.output().path
        )
        run_cmd(cmdline, dry_run=self.dry_run)

class vep_part(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return vt_part(infodict=self.infodict, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vt.vcf',
                                                           '.vep.vcf.gz'))

    def run(self):
        cmdline = """{vep} -i {vt_vcf} -o {vep_output_vcf} --vcf --cache --merged --fasta {REF} --sift b --polyphen b --symbol --numbers --biotype \
        --total_length --canonical --ccds --gene_phenotype --uniprot \
        --force_overwrite --offline --domains --regulatory --protein --tsl --variant_class --fork {threads} --force \
        --no_stats >> {vep_log} 2>&1""".format(
            vep=config.vep_pro,
            vt_vcf=self.input().path,
            REF=config.REF_file_path,
            vep_output_vcf=self.input().path.replace('.vt.vcf', '.vep.vcf'),
            threads=20,
            vep_log=self.input().path.replace('.vt.vcf', '.vep.log'))
        run_cmd(cmdline, dry_run=self.dry_run)

        cmdline = '{bgzip} -c {vep_output_vcf} > {vep_output_vcf_gz}'.format(
            bgzip=config.bgzip_pro,
            vep_output_vcf=self.input().path.replace('.vt.vcf', '.vep.vcf'),
            vep_output_vcf_gz=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)

        cmdline = '{tabix} -p vcf {vep_output_vcf_gz}'.format(
            tabix=config.tabix_pro,
            vep_output_vcf_gz=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class gemini_part(luigi.Task):
    infodict = luigi.DictParameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return vep_part(infodict=self.infodict,
                        dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vep.vcf.gz',
                                                           '.gemini.db'))

    def run(self):
        cmdline = """{gemini} load --cores {threads} -t VEP -v {vep_output_vcf_gz} {Output_db}; \
        {gemini} annotate -f {vep_output_vcf_gz} -a extract \
        -c SAD,SAF,AF,AD,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR \
        -t text,float,float,text,float,float,float,float,float \
        -o list,list,list,list,mean,mean,mean,mean,mean {Output_db} >> {gemini_log} 2>&1""".format(
            gemini=config.gemini_pro,
            threads=20,
            vep_output_vcf_gz=self.input().path,
            Output_db=self.output().path,
            gemini_log=self.output().path + '.log')
        run_cmd(cmdline, dry_run=self.dry_run)
