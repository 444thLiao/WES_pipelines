

"""
python ~/project/Whole_pipelines/aft_pipelines_analysis/filter_pipelines.py -pp 0,3,4,6 -o /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-57_all_filters_EXCEPT_AF_PASS_FILTER.csv

python ~/project/Whole_pipelines/aft_pipelines_analysis/filter_pipelines.py -pp 0,3,4,5,6 -o /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-57_all_filters_EXCEPT_AF_FILTER.csv

# upload
# filter coverage by hand.


python ~/project/Whole_pipelines/aft_pipelines_analysis/csv2bed.py /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-57_all_filters_EXCEPT_AF_PASS_FILTER.csv /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/bed/XK-57_all_filters_EXCEPT_AF_PASS_FILTER.bed

python ~/project/Whole_pipelines/aft_pipelines_analysis/extracted_pos_from_vcf.py /home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25.mt2.vcf.gz /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/bed/XK-25_filtered_somatic.bed /home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25T.mt2.vcf.gz /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vcf

vt decompose -s /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vcf.gz | vt normalize -r /home/liaoth/data/hg19/ucsc.hg19.fasta - > /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vt.vcf

mkdir /home/liaoth/project/180104_XK/gpz_server/gemini_annotate/pcgr/XK-25_FINAL

sudo python ~/data2/pcgr-0.5.3/pcgr.py --input_vcf /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vt.vcf  ~/data2/pcgr-0.5.3 '/home/liaoth/project/180104_XK/gpz_server/gemini_annotate/pcgr/XK-25_FINAL' ~/data2/pcgr-0.5.3/data/pcgr_configuration_default.toml XK-25_FINAL --force_overwrite
"""

