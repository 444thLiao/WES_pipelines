# Setting 结构

**`setting文件`理论上应该存在两个同样的拷贝，一个位于本地项目文件夹下，另一个位于服务器的项目文件夹下。**

<Br>总的来说分成三块

###1. 随着项目变更经常变更的参数
```python
##############################################################
##
## Usually you need to change part.
##############################################################

##Noramlly need to change part
###    Change Change
base_outpath = '/home/liaoth/project/XK_WES/180309_all/output'
base_inpath = '/home/liaoth/data_bank/XK_WES/*/'
bed_file_path = '/home/liaoth/data_bank/XK_WES/Sureselect_V6_COSMIC_formal.bed'
####For somatic pipelines
NORMAL_SIG = 'W'
TUMOR_SIG = 'T'
####Pair or not
Pair_data = True
###### programe path
vt_pro = '/home/liaoth/tools/vt/vt'
vep_pro = '/home/liaoth/tools/ensembl-vep-release-91/vep'
bgzip_pro = '/usr/bin/bgzip'
tabix_pro = '/usr/bin/tabix'
gemini_pro = '/usr/local/bin/gemini'
annovar_pro = '/home/liaoth/tools/annovar'
#####Paired format.
PE1_fmt = '{input}_R1'
PE2_fmt = '{input}_R2'
#####Single format
SE_fmt = '{input}_L001_R1_001'
######PCR gene number   ###PCR part
total_gen = 26
PCR_ON = False
###### luigi config
worker = 6
annovar_thread = 20
bip=''
max_memory= 10240
###### For parse file name in order to pipeline analyze.
sample_ID_pattern = 'XK-([0-9]+)[%s%s](?:-2)?' % (NORMAL_SIG, TUMOR_SIG)
pair_ID_pattern = '_R([12])'
mt2_for_pattern = 'XK-[0-9]+([%s%s])' % (NORMAL_SIG, TUMOR_SIG)
pair_name_pattern = '(XK-[0-9]+)[%s%s](?:-2)?' % (NORMAL_SIG, TUMOR_SIG)
sample_name_pattern = '(XK-[0-9]+[A-Z](?:-2)?)'


input_pair_tuple = [('XK-2W','XK-2T'),
                    ('XK-2W-2','XK-2T'),
                    ('XK-8W','XK-8T'),
                    ('XK-8W-2','XK-8T'),
                    ('XK-25W','XK-25T'),
                    ('XK-27W','XK-27T'),
                    ('XK-32W','XK-32T'),
                    ('XK-57W','XK-57T')]

PROJECT_NAME = 'XK'
```
可以通过留意`#`的注释，了解这一大块内部更小的结构。<Br>
需要注意的是，`base_inpath`的内容可以带有通配符(wildcard)，这是为了适应在**主流程**中的`auto_detect`的参数所进行的适应。

另外最后一个注释的信息，即各种正则表达式的模式(pattern)中，input_pair_tuple是为了更为准确的匹配复杂的情况，例如该例子中，存在`XK-2W-2 与XK-2W`两个正常组织的测序文件，并且同时与`XK-2T`成对，所以是比较复杂的情况，按照一般的识别过程可能无法正确的匹配，所以需要进行特殊的声明。**如果无复杂情况，可留空(即[])**。

###2.视情况而定的自适应参数

由于很多时候，给予的测序文件名难以找到一个共有的规律，并且有很多无用和容易混淆的字符，所以为了解决这个问题，我们添加了这个部分，通过简单的参数，来使得可以应对更加复杂的文件名

```python
##############################################################
#####
#####    SELF ADJUST PART OF INPUT FILENAME
##### USED IT IN SOME COMPLEXITY FILENAME YOU HARD TO PARSE
##############################################################
self_adjust_fn = True
filter_str = 'R-'
R1_INDICATOR = '_R1'
R2_INDICATOR = '_R2'
fq_suffix = '.fastq.gz'
```

正如前文所说，如果将`self_adjust_fn`设为`False`，那么后面的内容都不会生效。


###3. 固定的目录结构与数据库的位置

这一部分基本是不需要进行变动的，特别是目录结构。但是数据库可能会因为后续分析的要求的变化，包括测序方所使用的扩增试剂盒的不同或者要求的不同，有可能从hg19 -> hg38

```python

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

## DB files, Normaly don't need to change.
REF_file_path = '/home/db_public/hg19/ucsc.hg19.fasta'
known_gold_cvf = '/home/db_public/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
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
```
