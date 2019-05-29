# setting file , we normally change.

##############################################################
##
## Usually you need to change part.
##############################################################

##Noramlly need to change part
###    Change Change
base_outpath = '/home/liaoth/project/XK_WES/180321/output'
base_inpath = '/home/liaoth/data_bank/XK_WES/180321/'
bed_file_path = '/home/liaoth/data_bank/XK_WES/Sureselect_V6_COSMIC_formal.bed'
####For somatic pipelines
NORMAL_SIG = 'W'
TUMOR_SIG = 'T'
####Pair or not
Pair_data = True
###### server programe path
vt_pro = '/home/liaoth/tools/vt/vt'
vep_pro = '/home/liaoth/tools/ensembl-vep-release-91/vep'
bgzip_pro = '/usr/bin/bgzip'
tabix_pro = '/usr/bin/tabix'
gemini_pro = '/usr/local/bin/gemini'
annovar_pro = '/home/liaoth/tools/annovar'
gatk_pro = '/tools/GATK-4.0/gatk'
samtools_pro = '/home/liaoth/.local/bin/bin/samtools'
trimmomatic_jar = '/home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar'
pircard_jar = "/home/liaoth/tools/picard-tools-2.5.0/picard.jar"
gatkv36_path = "/home/liaoth/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"

server_script_path = '/home/liaoth/tools/Whole_pipelines'
###### Local setting.
server_path = 'liaoth@10.10.1.64'
local_project_path = '/home/liaoth/data2/project/XK_WES/180321/output'
local_vt_pro = '/home/liaoth/tools/vt/vt'
pcgr_pro = '/home/liaoth/data2/pcgr-0.5.3'
local_hg19_ref_db = '/home/liaoth/data2/hg19/ucsc.hg19.fasta'
local_bed_file = '/home/liaoth/data2/project/XK_WES/Sureselect_V6_COSMIC_formal.bed'
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

worker = 6

bip = ''
max_memory = 10240

## software params
trimmomatic_thread = 10
sort_sam_ram = "60G"
sort_sam_thread = 10
gatk_thread = 20
annovar_thread = 20
## DB files, Normaly don't need to change.
annovar_db = "/home/liaoth/tools/annovar/humandb/"
REF_file_path = '/home/db_public/hg19/ucsc.hg19.fasta'
known_gold_vcf = '/home/db_public/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
db_snp = '/home/db_public/hg19/dbsnp_138.hg19.vcf'
cos_snp = '/home/db_public/hg19/combinded_hg19_cosmic_v77.vcf'
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
