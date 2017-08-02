from collections import defaultdict

import luigi
from main import *
from special_fun import Add_cov_ino_in_vcf as P_vcf


class QC_trimmomatic(luigi.Task):
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)

    def output(self):
        sample_name = pfn(self.PE1, 'sample_name')
        project_name = pfn(self.PE1, 'project_name')
        return luigi.LocalTarget(
            '{base}/{PN}_result/trim_result/{SN}_trimed.log'.format(base=base_outpath, PN=project_name, SN=sample_name))

    def run(self):

        project_name = pfn(self.PE1, 'project_name')

        if os.path.isdir('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name)) != True:
            os.makedirs('{base}/{PN}_result/trim_result'.format(base=base_outpath, PN=project_name))

        input1 = self.PE1
        input2 = self.PE2
        if input2 != None:
            os.system(
                "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 {base_in}/{input1}.fastq.gz {base_in}/{input2}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz {base_out}/{input1}.unpaired.fq.gz {base_out}/{input2}.clean.fq.gz {base_out}/{input2}.unpaired.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
                    input1=input1, input2=input2, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                    output=self.output().path))
        else:
            os.system(
                "java -jar ~/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 {base_in}/{input1}.fastq.gz -trimlog {output} {base_out}/{input1}.clean.fq.gz ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                    input1=input1, base_in=base_inpath, base_out=self.output().path.rpartition('/')[0],
                    output=self.output().path))


class GenerateSam_pair(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        if Pair_data:
            input1 = PE1_fmt.format(input=self.sampleID)
            input2 = PE2_fmt.format(input=self.sampleID)
            return QC_trimmomatic(PE1=input1, PE2=input2)
        else:
            input1 = self.sampleID
            return QC_trimmomatic(PE1=input1)

    def output(self):
        sample_name = self.sampleID
        project_name = pfn(self.sampleID, 'project_name')

        return luigi.LocalTarget(
            output_fmt.format(path=base_outpath, PN=project_name, SN=sample_name) + '.sam')

    def run(self):
        sample_name = self.sampleID
        project_name = pfn(self.sampleID, 'project_name')

        if Pair_data:
            input1 = '{base_in}/{pe1_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe1_fmt=PE1_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            input2 = '{base_in}/{pe2_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                              pe2_fmt=PE2_fmt.format(input=self.sampleID),
                                                              input=self.sampleID)
            if not os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) == True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            os.system(
                "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1} {i2}  > {o}".format(
                    SN=sample_name, REF=REF_file_path, i1=input1, i2=input2, o=self.output().path))
        else:
            input1 = '{base_in}/{SE_fmt}.clean.fq.gz'.format(base_in=self.input().path.rpartition('/')[0],
                                                             SE_fmt=SE_fmt.format(input=self.sampleID),
                                                             input=self.sampleID)
            if os.path.isdir(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name)) != True:
                os.makedirs(output_dir.format(path=base_outpath, PN=project_name, SN=sample_name))
            os.system(
                "bwa mem -M -t 20 -k 19 -R '@RG\\tID:{SN}\\tSM:{SN}\\tPL:illumina\\tLB:lib1\\tPU:L001' {REF} {i1}  > {o}".format(
                    SN=sample_name, REF=REF_file_path, i1=input1, o=self.output().path))


class Convertbam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [GenerateSam_pair(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.sam','.bam'))

    def run(self):
        os.system("samtools view -F 0x100 -bSu %s -o %s" % (self.input()[0].path, self.output().path))


class sorted_bam(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [Convertbam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.bam','_sorted.bam'))

    def run(self):
        os.system("samtools sort -m 100G -f -@ 30 %s %s" % (self.input()[0].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


#########2
class MarkDuplicate(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('_sorted.bam','.dedup.bam'))

    def run(self):
        if PCR_ON:
            os.system('touch %s ' % self.output().path)
        else:
            os.system(
                "java -Xmx2g -jar ~/tools/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s/dedup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=true AS=true" % (
                    self.input()[0].path, self.output().path, self.output().path.rpartition('/')[0]))


#########3
class RealignerTargetCreator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID),sorted_bam(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam','.realign.intervals'))

    def run(self):
        if PCR_ON:
            os.system(
            "java -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R %s -I %s --known %s -o %s" % (
                REF_file_path, self.input()[0].path, known_gold_cvf, self.output().path))
        else:
            os.system(
            "java -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R %s -I %s --known %s -o %s" % (
                REF_file_path, self.input()[1].path, known_gold_cvf, self.output().path))


#########4
class IndelRealigner(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MarkDuplicate(sampleID=self.sampleID), RealignerTargetCreator(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.dedup.bam','.realign.bam'))

    def run(self):
        os.system(
            "java -Xmx5g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))
        os.system('samtools index %s' % self.output().path)


#########5
class BaseRecalibrator(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.rpartition('.realign.bam')[0] + '.recal_data.table')

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 25 -R %s -I %s -knownSites %s -knownSites %s -o %s" % (
                REF_file_path, self.input()[0].path, db_snp, known_gold_cvf, self.output().path))


#########6
class PrintReads(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [IndelRealigner(sampleID=self.sampleID), BaseRecalibrator(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.realign.bam','.recal_reads.bam'))

    def run(self):
        os.system(
            "java -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T PrintReads -R %s -I %s -BQSR %s -o %s" % (
                REF_file_path, self.input()[0].path, self.input()[1].path, self.output().path))
        os.system('samtools index %s' % self.output().path)



#########somatic pipeline
class MuTect2_pair(luigi.Task):
    ####combine N & T
    sample_IDs = luigi.Parameter()

    def requires(self):
        sampleIDs = self.sample_IDs.split(',')
        return [PrintReads(sampleID=sampleIDs[0]),PrintReads(sampleID=sampleIDs[1])]

    def output(self):
        sampleIDs = self.sample_IDs.split(',')
        Project_ID = pfn(sampleIDs[0], 'project_name')
        pair_name = pfn(sampleIDs[0], 'pair_name')

        output_path = somatic_pair_output_fmt.format(path=base_outpath,PN=Project_ID,PairN=pair_name)+'.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        sampleIDs = self.sample_IDs.split(',')

        output_dir = self.output().path.rpartition('/')[0]

        if os.path.isdir(output_dir) != True:
            os.makedirs(output_dir)

        if pfn(sampleIDs[0], 'mt2_for') == NORMAL_SIG:
            input_normal = self.input()[0].path
            input_tumor = self.input()[1].path
        elif pfn(sampleIDs[0], 'mt2_for') == TUMOR_SIG:
            input_normal = self.input()[1].path
            input_tumor = self.input()[0].path
        else:
            input_tumor=''
            input_normal=''

        prefix = self.output().path.rpartition('.bam')[0]

        os.system(
            '''java -Xmx10g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T MuTect2 --allSitePLs -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:normal {input_normal} --input_file:tumor {input_tumor} --out {prefix}.vcf --bamOutput {prefix}.bam --log_to_file {prefix}.log'''.format(
                REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_tumor=input_tumor, input_normal=input_normal,
                prefix=prefix))


class MuTect2_single(luigi.Task):
    sample_NT = luigi.Parameter()

    def requires(self):
        return [PrintReads(sampleID=self.sample_NT)]

    def output(self):

        Project_ID = pfn(self.sample_NT, 'project_name')
        sample_name = pfn(self.sample_NT, 'sample_name')

        output_path = somatic_single_output_fmt.format(path=base_outpath, PN=Project_ID, SN=sample_name) + '.mt2.bam'
        return luigi.LocalTarget(output_path)

    def run(self):
        input1 = self.input()[0].path
        mt2_id = pfn(self.sample_NT, 'mt2_for')
        prefix = self.output().path.rpartition('.bam')[0]
        output_dir = self.output().path.rpartition('/')[0]

        if os.path.isdir(output_dir) != True:
            os.makedirs(output_dir)

        if mt2_id == NORMAL_SIG:
            os.system(
                '''java -Xmx10g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T MuTect2 --allSitePLs --artifact_detection_mode -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_tumor} --out {prefix}.vcf --bamOutput {prefix}.bam --log_to_file {prefix}.log '''.format(
                    REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_tumor=input1, prefix=prefix))
        ####Normal only
        else:
            os.system(
                '''java -Xmx10g -jar ~/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T MuTect2 --allSitePLs --artifact_detection_mode -R {REF} --cosmic {cosmic} --dbsnp {db_snp} --input_file:tumor {input_tumor} --out {prefix}.vcf --bamOutput {prefix}.bam --log_to_file {prefix}.log --tumor_lod 4 '''.format(
                    REF=REF_file_path, cosmic=cos_snp, db_snp=db_snp, input_tumor=input1, prefix=prefix))
        ####Tumor only

class Add_cov_infos_SO(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [MuTect2_single(sample_NT=self.sampleID),
                PrintReads(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.mt2.bam', '.mt2.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_SO(self.input()[1].path, self.input()[0].path.replace('.bam','.vcf'), self.output().path, REF_file_path)


class Add_cov_infos_PA(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        _IDs = self.sampleID.split(',')
        if pfn(_IDs[0],'mt2_for') == NORMAL_SIG:
            normal_ID = _IDs[0]
            tumor_ID = _IDs[1]
        else:
            tumor_ID = _IDs[0]
            normal_ID = _IDs[1]

        return [MuTect2_pair(sample_IDs=self.sampleID),
                PrintReads(sampleID=normal_ID),
                PrintReads(sampleID=tumor_ID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.mt2.bam', '.mt2.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_PA([self.input()[1].path,self.input()[2].path],
                            self.input()[0].path.replace('.bam','.vcf'),
                            self.output().path,
                            REF_file_path)


#########15
class vt_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        if self.sampleID in pair_bucket:
            pair_value = pair_bucket[self.sampleID]
            return [Add_cov_infos_PA(sampleID=','.join(pair_value))]

        else:
            return [Add_cov_infos_SO(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.added_cov.vcf', '.vt.vcf'))

    def run(self):
        os.system("{vt} decompose -s {input_vcf} | {vt} normalize -r {REF} - > {vt_vcf}".format(
            vt=vt_pro, input_vcf=self.input()[0].path, REF=REF_file_path, vt_vcf=self.output().path))


class vep_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [vt_part(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.vt.vcf', '.vep.vcf.gz'))

    def run(self):
        os.system("""{vep} -i {vt_vcf} --cache --merged --fasta {REF} --sift b --polyphen b --symbol --numbers --biotype \
        --total_length --canonical --ccds -o {vep_output_vcf} --vcf --hgvs --gene_phenotype --uniprot \
        --force_overwrite --port 3337 --domains --regulatory --protein --tsl --variant_class --fork {threads} --force \
        --no_stats >> {vep_log} 2>&1""".format(vep=vep_pro,
                                               vt_vcf=self.input()[0].path,
                                               REF=REF_file_path,
                                               vep_output_vcf=self.input()[0].path.replace('.vt.vcf', '.vep.vcf'),
                                               threads=20,
                                               vep_log=self.input()[0].path.replace('.vt.vcf', '.vep.log')))
        os.system('{bgzip} -c {vep_output_vcf} > {vep_output_vcf_gz}'.format(bgzip=bgzip_pro,
                                                                             vep_output_vcf_gz=self.output().path,
                                                                             vep_output_vcf=self.input()[0].path.replace('.vt.vcf', '.vep.vcf')))
        os.system('{tabix} -p vcf {vep_output_vcf_gz}'.format(tabix=tabix_pro,
                                                              vep_output_vcf_gz=self.output().path))


class gemini_part(luigi.Task):
    sampleID = luigi.Parameter()

    def requires(self):
        return [vep_part(sampleID=self.sampleID)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.vep.vcf.gz', '.gemini.db'))

    def run(self):

        os.system("""{gemini} load --cores {threads} -t VEP -v {vep_output_vcf_gz} {Output_db}; {gemini} annotate -f \
        {vep_output_vcf_gz} -a extract -c SAD,SAF,AF,AD,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR -t text,text,text,text, \
        float,float,float,float,float -o list,list,list,list,mean,mean,mean,mean,mean {Output_db} >> {gemini_log} \
        2>&1""".format(
            gemini=gemini_pro,
            threads=20,
            vep_output_vcf_gz=self.input()[0].path,
            Output_db=self.output().path,
            gemini_log=self.output().path + '.log'))


class workflow(luigi.Task):
    x = luigi.Parameter()
    def requires(self):
        samples_IDs = str(self.x).split(',')

        pair_bucket = defaultdict(list)
        for _x in samples_IDs:
            pair_bucket[pfn(_x,'pair_name')].append(_x)
        global pair_bucket
        ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}

        samples_IDs += [_x for _x in pair_bucket.keys()]
        for i in samples_IDs:
            yield gemini_part(sampleID=i)




# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
