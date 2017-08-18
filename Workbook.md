##WES pipeliens 'server' & 'station'  workbook

####1. preprocessing
1. Figure out the format rules of fasta file.
2. Fullfill the [setting.py](/home/liaoth/project/Whole_pipelines/setting.py) file.*This is ori file. Maybe after updated field will not show in here.*Normally,you just need to **change these field**.
```python
base_outpath = '/home/liaoth/project_formal/170801_XK/output'
base_inpath = '/media/vd1/170801_X-TEN/C170010-P002/CHG024757/'
###### For parse file name in order to pipeline analyze.
sample_ID_pattern = 'XK-([0-9]{2})[%s%s]' % (NORMAL_SIG, TUMOR_SIG)
pair_ID_pattern = '_R([12])'
mt2_for_pattern = 'XK-[0-9]{2}([%s%s])' % (NORMAL_SIG, TUMOR_SIG)
pair_name_pattern = '(XK-[0-9]{2})[%s%s]' % (NORMAL_SIG, TUMOR_SIG)
sample_name_pattern = '(XK-[0-9]{2}[A-Z]-?[0-9]?)'
input_pair_tuple = []
PROJECT_NAME = 'XK'
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
3. Extract the sample ID which represents the data.e.g.*XK-8T-2,XK-8-2*
> Output a `setting.py`


- - -

####2. Run analysis pipelines and finish the server jobs.
#####2.1 run somatic pipelines
Using `python ~/tools/Whole_pipelines/main.py -A somatic -arg XK-25T,XK-25W,XK-27T,XK-27W -set /home/liaoth/project_formal/170801_XK/setting.py` to init the pipelines.**remember to disown the background jobs in order to hang up the task.**

#####2.2 run germline pipelines
Maybe you also need to `python ~/tools/Whole_pipelines/main.py -A germline -arg XK-25T,XK-25W,XK-27T,XK-27W -set /home/liaoth/project_formal/170801_XK/setting.py`. This is base on your need.

#####2.3 generate coverage info file
After `recal_reads.bam` is generated. You need to run `python ~/tools/Whole_pipelines/special_fun/cal_Cov_script_version.py -b {bam} -B  /home/liaoth/project_formal/170602_XK/WES_sureselect_v6.bed -r /home/liaoth/data/hg19/ucsc.hg19.fasta`. This is basically cal depth in all needed pos in WES project.
	**You also can use this shell script**
```shell
for path in `ls /home/liaoth/project_formal/170801_XK/output/XK_result/*/*.recal_reads.bam`; do python ~/tools/Whole_pipelines/special_fun/cal_Cov_script_version.py -b {bam} -B {bed} -r {ref} & done
for path in `ls /home/liaoth/project_formal/170801_XK/output/XK_result/*/*_sorted.bam`; do python ~/tools/Whole_pipelines/special_fun/cal_Cov_script_version.py -b {bam} -B {bed} -r {ref} & done
```
**After above process , you can get '.info' file.**

#####2.4 generate assessment table
1. Using **'.info' file ,raw fasta, trimed fasta** to generate **Accessment excel**.Example format see [accessment excel](/home/liaoth/project/170602_XK/WES两次活检总体评估表_改正.xlsx). Script is comming soon.
2. Using **.info file** to add coverage info into **pipelines generated csv**.Using `python ~/tools/Whole_pipelines/special_fun/add_per_info_into_csv.py -i ~/project_formal/170602_XK/output/XK_result/XK-8_somatic/ -o ~/project_formal/170602_XK/output/XK_result/XK-8_somatic/XK-8-2_with_info.csv -tb ~/project_formal/170602_XK/output/XK_result/XK-8T-2/XK-8T-2.recal_reads.bam -nb ~/project_formal/170602_XK/output/XK_result/XK-8W-2/XK-8W-2.recal_reads.bam`
**Of course you also can use bash script order plus this command to batch analysis**

#####2.5 Download all needed files to station
From station to download some file. I make a example table with one sample.e.g. XK-8
file suffix|file description|file location
----|----|----|
_with_info.csv|ALL snp csv plus corresponding coverage info|tumor single and pair,both.Total 2 file.
.merged.anno.hg19_multianno.csv|ALL snp CSV|above plus germline. Total is 3 file.
.vcf|ALL snp VCF file|Total is 4 file.
> Output a series of analysis file(.bam,.vcf,.bai,.csv,.info). 
> Output a assesment file.(**excel**)(Needed to **send first** in case some sample need to re-sequencing)

- - -

####3. Station workbook: after downloading all needed files.
#####3.1 Downloading all files
#####3.2 run filtering pipelines
Self edit [pipelines Script](/home/liaoth/project/Whole_pipelines/aft_pipelines_analysis/filter_pipelines.py) to run pipelines.
>Caution: You need to run twice each sample with different paramaters. One is without AF, another is without AF & PASS filter.
>Caution: After pipelines finished, it will output stats in this pipelines runing. you need to add it to a table.

**After this process, you will get two filtered csv file.
#####3.3 Drawing heatmap with coverage info.
You can use **2.3 generated info files** to draw the heatmap in specific genes.The gene list is in the [genes script](/home/liaoth/projcet/Whole_pipeliens/aft_pipelines_analysis/genes_list.py). Till 2017/08/15, We have four different set of genes. 
1. One is significance gene in lung cancer according to some paper. 
2. it is immune relative gene from Dong Hui. Three of them.
3. it is set of genes relative to dna repair mechanism.
4. clearSeq target gene. 151 gens.

Because of the different set of genes, we need to subset the big(~4G) coverage df into small one. You can use this [Script](/home/liaoth/project/Whole_pipelines/special_fun/cov_and_depths_heatmap.py) as references. But it can not use as auto script. you also need to self-adjust it in ipython.

#####3.4 Drawing venn graph with different sample different time.
Function have been written in [Script](/home/liaoth/project/Whole_pipelines/aft_pipelines_analysis/draw_venn.py). This Script Can help you draw 2-file compare venn and adjust the fontsize. This also can extract shared SNP presented in venn.


#####3.5 more Annotate snp with pcgr.
a. Extract snp info from csv as bed file.[function script](/home/liaoth/project/Whole_pipelines/aft_pipelines_analysis/csv2bed.py)
b. use bed to extract specific SNP from vcf.  [function script](/home/liaoth/project/Whole_pipelines/aft_pipelines_analysis/extracted_pos_from_vcf.py)
> You need to notice that, we have two csv file. So that is not easy to merge them together and extract specific pos. So we write this script to help guys do this jobs.
c. use pcgr to annote.
d. pack up infomation from pcgr.




