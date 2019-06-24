from __future__ import print_function
import glob, os, sys
import pandas as pd

"""
python 

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

script_path = os.path.abspath(__file__)
dir_script = os.path.dirname(script_path)

if len(sys.argv) == 3 and '/' in sys.argv[1]:
    setting_file = os.path.abspath(sys.argv[1])
    dir_path = os.path.dirname(setting_file)
    sys.path.insert(0, dir_path)
    from setting import *
else:

    setting_file = ''
    exit('Please using `aft_pipelines_main.py setting.py 1,4,5,6,7,8,9`.The setting.py should in your project path.')
    from setting import *

# relative input csv path
somatic_csv = os.path.join(local_project_path, 'somatic/')
germline_csv = os.path.join(local_project_path, 'germline/')
vcf_path = os.path.join(local_project_path, 'vcf_storge/')
info_summary_path = os.path.join(local_project_path, 'summary_info/')
# relative output path
final_csv = os.path.join(local_project_path, 'final_output/')
filtered_csv = os.path.join(local_project_path, 'filtered_file/')
pcgr_output = os.path.join(local_project_path, 'pcgr/')
summary_depth = os.path.join(local_project_path, 'info_summary/')
pack_result = os.path.join(local_project_path, '../packet_result')


def check_path(path):
    if not os.path.isdir(path):
        os.system('mkdir -p %s ' % path)


for opath in [somatic_csv, germline_csv, vcf_path, final_csv, filtered_csv, pcgr_output, pack_result,info_summary_path]:
    check_path(opath)


def run(cmd):
    print(cmd)
    os.system(cmd)


def download_scp(only_info_summary=False):
    if only_info_summary:
        # target cov_info_summary
        cmdline = 'scp %s:%s/%s_result/*/*_cov_summary.info %s' % (
            server_path, base_outpath, PROJECT_NAME, info_summary_path)
        run(cmdline)
        exit('Finished download all need file.')
    # target csv file
    cmdline = 'scp %s:%s/%s_result/*_somatic/*anno.csv %s' % (server_path, base_outpath, PROJECT_NAME,
                                                              somatic_csv)
    run(cmdline)
    cmdline = 'scp %s:%s/%s_result/*W/*anno.csv %s' % (server_path, base_outpath, PROJECT_NAME, germline_csv)
    run(cmdline)
    # target vcf files
    cmdline = "scp %s:%s/%s_result/*_somatic/*.mt2.vcf %s" % (server_path, base_outpath, PROJECT_NAME, vcf_path)
    run(cmdline)
    cmdline = 'scp %s:%s/%s_result/*W/*merged.vcf %s' % (server_path, base_outpath, PROJECT_NAME, vcf_path)
    run(cmdline)
    print('finished download all need file.')


def run_pipelines():
    from filter_pipelines import filter_pipelines2
    samples_ids = [os.path.basename(i).split('.mt2.merged.')[0] for i in glob.glob(os.path.join(somatic_csv, '*.csv'))]
    pair_name = [i for i in samples_ids if NORMAL_SIG not in i and TUMOR_SIG not in i]
    for each_pair in pair_name:
        if each_pair.count('-') == 2:
            normal_name = each_pair.rpartition('-')[0] + NORMAL_SIG
            tumor_name = each_pair.rpartition('-')[0] + TUMOR_SIG + '-' + each_pair.rpartition('-')[2]
        else:
            normal_name = each_pair + NORMAL_SIG
            tumor_name = each_pair + TUMOR_SIG
        germline = glob.glob(os.path.join(germline_csv, normal_name + '.merged*.csv'))[0]
        somatic_normal = glob.glob(os.path.join(somatic_csv, normal_name + '.mt2*.csv'))[0]
        somatic_tumor = glob.glob(os.path.join(somatic_csv, tumor_name + '.mt2*.csv'))[0]
        somatic_pair = glob.glob(os.path.join(somatic_csv, each_pair + '.mt2*.csv'))[0]
        output_file = os.path.join(filtered_csv, '%s_all_except_AF_depth.csv' % each_pair)
        if not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        print("filter_pipelines2(%s,%s,%s,%s,%s,pp=[3,4,5,6])" % (
            germline, somatic_normal, somatic_tumor, somatic_pair, output_file))
        filter_pipelines2(germline, somatic_normal, somatic_tumor, somatic_pair, output_file, pp=[3, 4, 5, 6])
        filter_pipelines2(germline, somatic_normal, somatic_tumor, somatic_pair,
                          output_file.replace('except_AF_depth.csv', 'except_AF_depth_PASS.csv'), pp=[3, 4, 6])

    print('finished filterered all samples')


def download_filtered_files(fetch_path):
    cmdline = "scp %s:%s %s" % (server_path, fetch_path, final_csv)
    run(cmdline)
    print('finished download all final csv.')


def filter_out_10(outpath, threshold=10):
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
    file_df.loc[high_cov_bucket, :].to_csv(outpath, index=False)


def pack_its_up():
    cmdline1 = 'mkdir -p %s/variants_csv; cp %s/*_all_except_AF_PASS_with_info.csv %s/variants_csv/' % (
        pack_result, final_csv, pack_result)
    cmdline2 = 'mkdir -p %s/html_report; cp %s/*/*.html %s/html_report/' % (pack_result, pcgr_output, pack_result)
    from simpy_summary import make_summary
    function_csv = make_summary(glob.glob('%s/variants_csv/*.csv' % pack_result))
    function_csv.to_csv('%s/variation_function_stats.csv' % pack_result, index=False)
    # cmdline3 = 'python3 %s/draw_quality_line.py %s' % (dir_script, project_name)
    run(cmdline1)
    run(cmdline2)
    # run(cmdline3)


def remind_text(local_project_path):
    os.chdir(os.path.dirname(local_project_path))
    exec ('from setting import *')

    server_setting_path = os.path.join(os.path.dirname(base_outpath.rstrip('/')), 'setting.py')
    remind_run_command = 'Command you may need to run at server. Listed below:\n'
    remind_run_command += '''
    ##run script in order to generate accessment file. \n\n \
    /home/liaoth/tools/pysamstats_venv/bin/python2.7 %s/pre_pipelines_analysis/quality_accessment.py %s \n''' % (server_script_path, server_setting_path)
    remind_run_command += '#' * 50 + '\n'
    remind_run_command += 'run cal_cov script to get coverage info from different bam files.'
    remind_run_command += '''for each in %s/XK_result/*/*sorted.bam; do python %s/pre_pipelines_analysis/cal_Cov_script_version.py -b $each -B %s -r %s & done''' % (
    base_outpath, server_script_path, bed_file_path, REF_file_path)
    remind_run_command += '''for each in %s/XK_result/*/*recal_reads.bam; do python %s/pre_pipelines_analysis/cal_Cov_script.py -b $each -B %s -r %s & done''' % (
    base_outpath, server_script_path, bed_file_path, REF_file_path)
    remind_run_command += '''
        ##run script in order to generate accessment file. \n\n \
        python %s/aft_pipelines_analysis/run_info_summary.py %s \n''' % (
    server_script_path, server_setting_path)
    remind_run_command += '''
    ##run script which is fetch cov_info from .info file and add it into csvfile. \n\n \
    python2 Whole_pipelines/aft_pipelines_analysis/run_add_per_info_into_csv.py %s \n''' % server_setting_path
    return remind_run_command


def draw_coverage_depths():
    from draw_quality_line import draw_coverage_depths as dcd
    dcd(info_summary_path, NORMAL_SIG, TUMOR_SIG)
    print('finish drawing.')


if __name__ == '__main__':
    args = sys.argv[-1].split(',')
    print('Using 1,4,5,6,7,8,9 as parameter to init all process. Or it will not do anything.')
    print('Downloading all files and processing all files in local path: %s' % (
        local_project_path))
    print('#' * 40)
    print('#' * 40)
    print(remind_text(local_project_path))
    print('#' * 40)
    print('#' * 40)

    server_setting_path = os.path.join(local_project_path, 'setting.py')
    if '1' in args:
        if str(input('Download ? Y/y')).upper() == 'Y':
            download_scp()
        else:
            print('wrong command,just pass.')
    if '2' in args:
        download_scp(only_info_summary=True)
    if '3' in args:
        draw_coverage_depths()
    if '4' in args:
        print('Runing filtered pipelines. Output all samples filtered file into %s' % filtered_csv)
        run_pipelines()
    if '5' in args:

        t_path = '%s/temp_/*_with_info.csv' % os.path.dirname(base_outpath.rstrip('/'))
        print(
            '''You need to upload these file to server to cal coverage based on bam files and run comands like this.\n\n Recommanded upload path: \nscp {filtered_csv}/*_all_except_AF_depth_PASS.csv {server_path}:{project_path}/temp_/. \nOr you will need to modify script 'run_add_per_info_into_csv.py'.\nAfter you finished, please pass the path with regular expression files to this.'''.format(
                project_path=os.path.dirname(base_outpath.rstrip('/')), filtered_csv=filtered_csv,
                server_path=server_path))
        if str(input('Make sure your path is %s. Y/y' % t_path)).upper() == 'Y':
            download_filtered_files(t_path)
        else:
            print(
                "Wrong input. Maybe you need to stop. If is 'underfind' error, you may need to use python3 utill you figure it out. Or you can directly ask author.")
    if '6' in args:
        print("*" * 50)
        print("Turning final csv file into bed.")
        print("*" * 50)
        for with_info_csv in glob.glob('%s/*_with_info.csv' % final_csv):
            filter_out_10(with_info_csv)
            cmdline = 'python %s/csv2bed.py %s %s' % (
                dir_script, with_info_csv, with_info_csv.replace('.csv', '.bed'))
            run(cmdline)
    if '7' in args:
        print("*" * 50)
        print("Using bed file to fetch pos from vcf file.")
        print("*" * 50)
        for bed_file in glob.glob('%s/*.bed' % final_csv):
            each_pair = os.path.basename(bed_file).split('_all')[0]
            if each_pair.count('-') == 2:
                normal_name = each_pair.rpartition('-')[0] + NORMAL_SIG
                tumor_name = each_pair.rpartition('-')[0] + TUMOR_SIG + '-' + each_pair.rpartition('-')[2]
            else:
                normal_name = each_pair + NORMAL_SIG
                tumor_name = each_pair + TUMOR_SIG
            cmdline = 'python %s/extracted_pos_from_vcf.py %s %s %s %s %s' % (
                dir_script,
                server_setting_path,
                bed_file,
                '%s/%s.mt2.vcf' % (vcf_path, each_pair),
                '%s/%s.mt2.vcf' % (vcf_path, tumor_name),
                '%s/%s_TS_merged.mt2.vcf' % (vcf_path, each_pair))
            run(cmdline)
    if '8' in args:
        print("*" * 50)
        print("Using fetched vcf to gernerate report base on pcgr.")
        print("*" * 50)

        for merged_vcf in glob.glob('%s/*_TS_merged.mt2.vcf.gz' % vcf_path):
            cmdline = '{vt_path} decompose -s {ori_vcf} | {vt_path} normalize -r {hg19_ref_path} - > {aft_vcf}'.format(
                ori_vcf=merged_vcf, aft_vcf=merged_vcf.replace('.vcf.gz', '.vt.vcf'),
                vt_path=local_vt_pro, hg19_ref_path=local_hg19_ref_db)
            run(cmdline)

        for vt_vcf in glob.glob('%s/*_TS_merged.mt2.vt.vcf' % vcf_path):
            pair_name = os.path.basename(vt_vcf).split('_TS')[0]
            cmdline = 'mkdir -p %s/%s_FINAL;' % (pcgr_output, pair_name)
            cmdline += "sudo python {pcgr_dir}/pcgr.py --input_vcf {in_vcf}  {pcgr_dir} '{pcgr_output}/{pair_name}_FINAL' {pcgr_dir}/data/pcgr_configuration_default.toml {pair_name}_FINAL --force_overwrite"
            cmdline = cmdline.format(in_vcf=vt_vcf, pcgr_output=pcgr_output, pair_name=pair_name, pcgr_dir=pcgr_pro)
            run(cmdline)
    print("*" * 50)
    print("Finishing all aft pipelines work. Congratulation~~~For Fun!")
    print("*" * 50)
    if '9' in args:
        pack_its_up()
