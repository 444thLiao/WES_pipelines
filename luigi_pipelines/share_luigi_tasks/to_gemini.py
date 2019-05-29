
import luigi


from special_fun import Add_cov_ino_in_vcf as P_vcf
from ... import run_cmd
from main import *

#########14
class Add_cov_infos(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        "abstract for implement"
        pass
        # return [CombineVariants(sampleID=self.sampleID, dry_run=self.dry_run),
        #         HaplotypeCaller(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        "abstract for implement"
        # return luigi.LocalTarget(self.input()[0].path.replace('.merged.vcf',
        #                                                       '.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_SO(self.input()[1].path,
                            self.input()[0].path,
                            self.output().path,
                            REF_file_path)


#########15
class vt_part(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return Add_cov_infos(sampleID=self.sampleID, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.added_cov.vcf',
                                                              '.vt.vcf'))

    def run(self):
        cmdline = "{vt} decompose -s {input_vcf} | {vt} normalize -r {REF}- > {vt_vcf}".format(
            vt=vt_pro,
            input_vcf=self.input().path,
            REF=REF_file_path,
            vt_vcf=self.output().path
        )
        run_cmd(cmdline, dry_run=self.dry_run)


class vep_part(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return vt_part(sampleID=self.sampleID, dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vt.vcf',
                                                           '.vep.vcf.gz'))

    def run(self):
        cmdline = """{vep} -i {vt_vcf} --cache --merged --fasta {REF} --sift b --polyphen b --symbol --numbers --biotype
        --total_length --canonical --ccds -o {vep_output_vcf} --vcf --hgvs --gene_phenotype --uniprot
        --force_overwrite --port 3337 --domains --regulatory --protein --tsl --variant_class --fork {threads} --force
        --no_stats >> {vep_log} 2>&1""".format(
            vep=vep_pro,
            vt_vcf=self.input().path,
            REF=REF_file_path,
            vep_output_vcf=self.input().path.replace('.vt.vcf', '.vep.vcf'),
            threads=20,
            vep_log=self.input().path.replace('.vt.vcf', '.vep.log'))
        run_cmd(cmdline, dry_run=self.dry_run)

        cmdline = '{bgzip} -c {vep_output_vcf} > {vep_output_vcf_gz}'.format(
            bgzip=bgzip_pro,
            vep_output_vcf=self.input().path.replace('.vt.vcf', '.vep.vcf'),
            vep_output_vcf_gz=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)

        cmdline = '{tabix} -p vcf {vep_output_vcf_gz}'.format(
            tabix=tabix_pro,
            vep_output_vcf_gz=self.output().path)
        run_cmd(cmdline, dry_run=self.dry_run)


class gemini_part(luigi.Task):
    sampleID = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)

    def requires(self):
        return vep_part(sampleID=self.sampleID,
                        dry_run=self.dry_run)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace('.vep.vcf.gz',
                                                           '.gemini.db'))

    def run(self):
        cmdline = """{gemini} load --cores {threads} -t VEP -v {vep_output_vcf_gz} {Output_db};
        {gemini} annotate -f {vep_output_vcf_gz} -a extract
        -c SAD,SAF,AF,AD,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR
        -t text,float,float,text,float,float,float,float,float
        -o list,list,list,list,mean,mean,mean,mean,mean {Output_db} >> {gemini_log} 2>&1""".format(
            gemini=gemini_pro,
            threads=20,
            vep_output_vcf_gz=self.input().path,
            Output_db=self.output().path,
            gemini_log=self.output().path + '.log')
        run_cmd(cmdline, dry_run=self.dry_run)