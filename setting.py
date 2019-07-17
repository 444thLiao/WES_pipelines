# setting file , we normally just verify the path of executor.

##############################################################
##
## Usually you need to change part.
##############################################################

##Noramlly need to change part
###    Change Change
from os.path import dirname, join
project_root_path = "/home/liaoth/project/Whole_pipelines/"

###### server programe path
bcftools_path = '/home/liaoth/tools/bcftools/bcftools'
vt_pro = '/home/liaoth/tools/vt/vt'
vep_pro = '/home/liaoth/tools/ensembl-vep/vep'
bgzip_pro = '/usr/bin/bgzip'
tabix_pro = '/usr/bin/tabix'
gemini_pro = '/usr/local/bin/gemini' # (optional, could stay default)
annovar_pro = '/home/liaoth/tools/annovar'
gatk_pro = '/home/liaoth/tools/gatk4/gatk'
samtools_pro = '/usr/bin/samtools'
samtools_version = 1.7  # difference between some versions
trimmomatic_jar = '/home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar'
pircard_jar = "/home/liaoth/tools/picard-tools-2.5.0/picard.jar"
gatkv36_path = "/home/liaoth/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
pcgr_dir = "/home/liaoth/tools/pcgr"
pcgr_toml_file = "/home/liaoth/tools/pcgr/conf/GPZ_pcgr_no_control.toml"
bed_file_path = join(project_root_path,
                     "db",
                     "Sureselect_V6_COSMIC_formal.bed")
######PCR gene number   ###PCR part
total_gen = 26
PCR_ON = False


##############################################################
##
## Usually you don't need to change part.
##############################################################


## file structure of output result, Normaly don't need to change.
output_dir = '{path}/{PN}_result/{SN}'
output_fmt = '{path}/{PN}_result/{SN}/{SN}'
somatic_pair_output_fmt = '{path}/{PN}_result/{PairN}_somatic/{PairN}'
somatic_single_output_fmt = '{path}/{PN}_result/{SN}_somatic/{SN}'
trim_fmt = '{base}/{PN}_result/trim_result'
filtered_dir = "{path}/filtered_csv/{PairN}"

bip = ''
max_memory = 10240

## software params
trimmomatic_thread = 5
sort_sam_ram = "2G"
sort_sam_thread = 2
gatk_thread = 5
annovar_thread = 1  # be carefull this... each one will take a lot of memory
gemini_thread = 5
vep_thread = 5
java_option = "-Xmx4g"
bwa_thread = 5

## DB files, Normaly don't need to change.
annovar_db = "/home/liaoth/tools/annovar/humandb/"
REF_file_path = '/home/liaoth/data/hg19/ucsc.hg19.fasta'
known_gold_vcf = '/home/liaoth/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
db_snp = '/home/liaoth/data/hg19/dbsnp_138.hg19.vcf'
cos_snp = '/home/liaoth/data/hg19/combinded_hg19_cosmic_v77.vcf'
genome_version = 'hg19'
db_names = 'refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,exac03,snp138,ljb26_all,clinvar_20160302'
####For cal_fun
cal_sample_name = ''
####For analysis after somatic pipelines.
snp138_common_file = '/home/liaoth/data/humandb/snp138Common.name.txt'
###### For format file. Normaly don't need to change.
columns_need_to_extract = ['A',
                           'C',
                           'G',
                           'T',
                           'A_Rate',
                           'C_Rate',
                           'G_Rate',
                           'T_Rate']

column_retained = ['Chr',
                   'Start',
                   'End',
                   'Ref',
                   'Alt',
                   'Func.refGene',
                   'Gene.refGene',
                   'GeneDetail.refGene',
                   'ExonicFunc.refGene',
                   'AAChange.refGene',
                   'A',
                   'C',
                   'G',
                   'T',
                   'A_Rate',
                   'C_Rate',
                   'G_Rate',
                   'T_Rate',
                   'snp138',
                   'CLINSIG',
                   'CLNDBN',
                   'CLNACC',
                   'CLNDSDB',
                   'CLNDSDBID',
                   'GT', 'AD', 'AF', 'Filter']

blastn_header = ['qseid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                 'evalue', 'bitscore']

bed_file_fmt = 'chr\tstart\tend\tgene_name\tstrand\n'
