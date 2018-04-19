from __future__ import print_function
import glob, os
import pandas as pd
"""
python ~/project/Whole_pipelines/aft_pipelines_analysis/filter_pipelines.py -i /home/liaoth/data2/project/XK_WES/180321/output/ -o /home/liaoth/data2/project/XK_WES/180321/output/final_output/ -s /home/liaoth/data2/project/XK_WES/180321/setting.py

# upload
# filter coverage by hand.


python ~/project/Whole_pipelines/aft_pipelines_analysis/csv2bed.py /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-57_all_filters_EXCEPT_AF_PASS_FILTER.csv /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/bed/XK-57_all_filters_EXCEPT_AF_PASS_FILTER.bed

python ~/project/Whole_pipelines/aft_pipelines_analysis/extracted_pos_from_vcf.py /home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25.mt2.vcf.gz /home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/bed/XK-25_filtered_somatic.bed /home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25T.mt2.vcf.gz /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vcf

vt decompose -s /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vcf.gz | vt normalize -r /home/liaoth/data/hg19/ucsc.hg19.fasta - > /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vt.vcf

mkdir /home/liaoth/project/180104_XK/gpz_server/gemini_annotate/pcgr/XK-25_FINAL

sudo python ~/data2/pcgr-0.5.3/pcgr.py --input_vcf /home/liaoth/data2/project/180104_XK/gpz_server/vcf_storge/final_vcf/XK-25_TS_merged.mt2.vt.vcf  ~/data2/pcgr-0.5.3 '/home/liaoth/project/180104_XK/gpz_server/gemini_annotate/pcgr/XK-25_FINAL' ~/data2/pcgr-0.5.3/data/pcgr_configuration_default.toml XK-25_FINAL --force_overwrite
"""
############################################################
## Local config plus some default config
############################################################
## 18/03/23
## absoulte path
project_name = '180321'
base_inpath = '/home/liaoth/data2/project/XK_WES/%s/output' % project_name

# relative input csv path
somatic_csv = os.path.join(base_inpath, 'somatic/')
germline_csv = os.path.join(base_inpath, 'germline/')
vcf_path = os.path.join(base_inpath, 'vcf_storge/')
# relative output path
final_csv = os.path.join(base_inpath, 'final_output/')
filtered_csv = os.path.join(base_inpath, 'filtered_file/')
pcgr_output = os.path.join(base_inpath,  'pcgr/')
summary_depth = os.path.join(base_inpath,'info_summary/')
pack_result = os.path.join(base_inpath,'../packet_result')

def check_path(path):
    if not os.path.isdir(path):
        os.system('mkdir -p %s ' % path)

for each in [somatic_csv,germline_csv,vcf_path,final_csv,filtered_csv,pcgr_output,pack_result]:
    check_path(each)
def run(cmd):
    print(cmd)
    os.system(cmd)


def download_scp(project_name):
    # target csv file
    cmdline = 'scp $lth_gpz_server:~/project/XK_WES/%s/output/XK_result/*_somatic/*anno.csv %s' % (project_name, somatic_csv)
    run(cmdline)
    cmdline = 'scp $lth_gpz_server:~/project/XK_WES/%s/output/XK_result/*W/*anno.csv %s' % (project_name, germline_csv)
    run(cmdline)
    # target vcf files
    cmdline = 'scp $lth_gpz_server:~/project/XK_WES/%s/output/XK_result/*_somatic/*.mt2.vcf %s' % (
    project_name, vcf_path)
    run(cmdline)
    cmdline = 'scp $lth_gpz_server:~/project/XK_WES/%s/output/XK_result/*W/*merged.vcf %s' % (project_name, vcf_path)
    run(cmdline)
    print('finished download all need file.')


def run_pipelines():
    cmdline = "python ~/project/Whole_pipelines/aft_pipelines_analysis/filter_pipelines.py -i %s -s %s -o %s" % (
    base_inpath,
    os.path.join(os.path.dirname(base_inpath), 'setting.py'), filtered_csv)
    run(cmdline)
    print('finished filterered all samples')


def download_filtered_files(fetch_path):
    cmdline = "scp $lth_gpz_server:%s %s" % (fetch_path,final_csv)
    run(cmdline)
    print('finished download all final csv.')

def filter_out_10(outpath,threshold=10):
    file_df = pd.read_csv(outpath)
    if 'N_mut_cov' not in file_df:
        return 'Wrong Columns'
    file_df.index = range(file_df.shape[0])
    file_df = file_df.loc[file_df.N_mut_cov.astype(str) != 'Off target', :]
    file_df = file_df.loc[file_df.T_mut_cov.astype(str) != 'Off target', :]
    _N_cov = file_df.loc[:, ['N_ref_cov', 'N_mut_cov']].astype(float).sum(1)
    _T_cov = file_df.loc[:, ['T_ref_cov', 'T_mut_cov']].astype(float).sum(1)
    high_cov_bucket = list(file_df.loc[(_N_cov >= threshold) & (_T_cov >= threshold) & (
            file_df.loc[:, 'T_mut_per'] > file_df.loc[:, 'N_mut_per']), :].index)
    file_df.loc[high_cov_bucket,:].to_csv(outpath,index=False)

def pack_its_up():
    cmdline1 = 'mkdir -p %s/variants_csv; cp %s/*_all_except_AF_PASS_with_info.csv %s/variants_csv/' % (pack_result,final_csv,pack_result)
    cmdline2 = 'mkdir -p %s/html_report; cp %s/*/*.html %s/html_report/' % (pack_result,pcgr_output,pack_result)
    from simpy_summary import make_summary
    function_csv = make_summary(glob.glob('%s/variants_csv/*.csv' % pack_result))
    function_csv.to_csv('%s/variation_function_stats.csv' % pack_result,index=False)
    cmdline3 = 'python3 /home/liaoth/project/Whole_pipelines/aft_pipelines_analysis/draw_quality_line.py %s' % project_name
    run(cmdline1)
    run(cmdline2)
    run(cmdline3)

def remind_text(base_inpath):
    #setting_py = os.path.join(os.path.dirname(base_inpath), 'setting.py')
    os.chdir(os.path.dirname(base_inpath))
    from setting import *
    server_setting_path = os.path.join(base_inpath,'setting.py')
    remind_run_command = ''
    remind_run_command += '''##/home/liaoth/tools/pysamstats_venv/bin/python2.7 /home/liaoth/tools/Whole_pipelines/pre_pipelines_analysis/quality_accessment.py %s \n''' % server_setting_path

if __name__ == '__main__':
    import sys
    args = sys.argv[-1].split(',')
    print('Using 1,2,3,4,5,6 as parameter to init all process. Or it will not do anything.')
    if '1' in args:
        if input('Download ? Y/y').upper() == 'Y':
            download_scp(project_name)
        else:
            print('wrong command,just pass.')

    if '2' in args:
        print('Runing filtered pipelines. Output all samples filtered file into %s' % filtered_csv)
        run_pipelines()
    if '3' in args:
        fetch_path = input('''You need to upload these file to server to cal coverage based on bam files.After you finished, please pass the path with regular expression files to this.''')
        if input('Make sure your path is %s . Y/y' % fetch_path).upper() == 'Y':
            download_filtered_files(fetch_path)
        else:
            print("Wrong input. Maybe you need to stop.Doesn't downlaod any things")
    if '4' in args:
        print("*" * 50)
        print("Turning final csv file into bed.")
        print("*" * 50)
        for each in glob.glob('%s/*_with_info.csv' % final_csv):
            filter_out_10(each)
            cmdline = 'python ~/project/Whole_pipelines/aft_pipelines_analysis/csv2bed.py %s %s' % (
            each, each.replace('.csv', '.bed'))
            run(cmdline)
    if '5' in args:
        print("*" * 50)
        print("Using bed file to fetch pos from vcf file.")
        print("*" *50)
        for each in glob.glob('%s/*.bed' % final_csv):
            each_pair = os.path.basename(each).split('_all')[0]
            if each_pair.count('-') == 2:
                normal_name = each_pair.rpartition('-')[0] + 'W'
                tumor_name = each_pair.rpartition('-')[0] + 'T' + '-' + each_pair.rpartition('-')[2]
            else:
                normal_name = each_pair + 'W'
                tumor_name = each_pair + 'T'
            cmdline = 'python ~/project/Whole_pipelines/aft_pipelines_analysis/extracted_pos_from_vcf.py %s %s %s %s' % (
            '%s/%s.mt2.vcf' % (vcf_path,each_pair), each,
            '%s/%s.mt2.vcf' % (vcf_path,tumor_name),
            '%s/%s_TS_merged.mt2.vcf' % (vcf_path,each_pair))
            run(cmdline)
    if '6' in args:
        print("*" * 50)
        print("Using fetched vcf to gernerate report base on pcgr.")
        print("*" *50)

        for each in glob.glob('%s/*_TS_merged.mt2.vcf.gz' % vcf_path):
            cmdline = '/home/liaoth/tools/vt/vt decompose -s %s | /home/liaoth/tools/vt/vt normalize -r /home/liaoth/data/hg19/ucsc.hg19.fasta - > %s' % (
            each, each.replace('.vcf.gz', '.vt.vcf'))
            run(cmdline)

        for each in glob.glob('%s/*_TS_merged.mt2.vt.vcf' % vcf_path):
            pair_name = os.path.basename(each).split('_TS')[0]
            cmdline = 'mkdir -p %s/%s_FINAL;' % (pcgr_output,pair_name)
            cmdline += "sudo python ~/data2/pcgr-0.5.3/pcgr.py --input_vcf %s  ~/data2/pcgr-0.5.3 '%s/%s_FINAL' ~/data2/pcgr-0.5.3/data/pcgr_configuration_default.toml %s_FINAL --force_overwrite"
            cmdline = cmdline % (each,pcgr_output, pair_name, pair_name)
            run(cmdline)
    print("*" * 50)
    print("Finishing all aft pipelines work. Congratulation~~~For Fun!")
    print("*" *50)
    if '7' in args:
        pack_its_up()
