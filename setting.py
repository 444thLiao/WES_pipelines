# setting file , we normally change.

##############################################################
##
## Usually you need to change part.
##############################################################

##Noramlly need to change part
###    Change Change
base_outpath = '/home/liaoth/project_formal/170123_XK/output'
base_inpath = '/home/shared/170123_new/fqs'
bed_file_path = ''
####For somatic pipelines
NORMAL_SIG = 'W'
TUMOR_SIG = 'T'
####Pair or not
Pair_data = True
###### programe path
vt_pro = '/home/liaoth/tools/vt/vt'
vep_pro = '/home/liaoth/tools/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl'
bgzip_pro = '/usr/bin/bgzip'
tabix_pro = '/usr/bin/tabix'
gemini_pro = '/usr/local/bin/gemini'
annovar_pro = '/home/liaoth/tools/annovar'
samtools_pro = '/home/liaoth/.local/bin/bin/samtools'
gatk_pro = '/tools/GATK-4.0/gatk'
#####Paired format.
PE1_fmt = '{input}_R1_001'
PE2_fmt = '{input}_R2_001'
#####Single format
SE_fmt = '{input}_R1_001'
######PCR gene number   ###PCR part
total_gen = 26
PCR_ON = True
###### For parse file name in order to pipeline analyze.
sample_ID_pattern = 'S([0-9]+)'
pair_ID_pattern = 'R([12])'
mt2_for_pattern = '[0-9]([%s%s])_S.*' % (NORMAL_SIG, TUMOR_SIG)
pair_name_pattern = '(.*)[%s%s]_S' % (NORMAL_SIG, TUMOR_SIG)
sample_name_pattern = '.*'


input_pair_tuple = ()
PROJECT_NAME = None


##############################################################
#####
#####    SELF ADJUST PART OF INPUT FILENAME
##### USED IT IN SOME COMPLEXITY FILENAME YOU HARD TO PARSE
##############################################################
self_adjust_fn = False
filter_str = ''
R1_INDICATOR = '_R1'
R2_INDICATOR = '_R2'
fq_suffix = '.fastq.gz'

##############################################################
##
## Usually you don't need to change part.
##############################################################

## file structure of output result, Normaly don't need to change.
output_dir = '{path}/{PN}_result/{SN}'
output_fmt = '{path}/{PN}_result/{SN}/{SN}'
somatic_pair_output_fmt = '{path}/{PN}_result/{PairN}_somatic/{PairN}'
somatic_single_output_fmt = '{path}/{PN}_result/{SN}_somatic/{SN}'


## For testing the protype of multithread
somatic_pair_output_fmt2 = '{path}/{PN}_result/{PairN}_somatic/{PairN}'
somatic_single_output_fmt2 = '{path}/{PN}_result/{SN}_somatic/{SN}'

worker = 12
annovar_thread = 20
bip='172.20.114.242'
max_memory= 10240
# unit is kb. default is 10G

## DB files, Normaly don't need to change.
REF_file_path = '/home/liaoth/data/hg19/ucsc.hg19.fasta'
known_gold_cvf = '~/data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
db_snp = '~/data/hg19/dbsnp_138.hg19.vcf'
cos_snp = '~/data/combinded_hg19_cosmic_v77.vcf'
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
