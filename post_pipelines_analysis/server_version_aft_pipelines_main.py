from __future__ import print_function
import glob, os, sys
import pandas as pd

############################################################
## Local config plus some default config
############################################################
## 18/04/24
## absoulte path

script_path = os.path.abspath(__file__)
dir_script = os.path.dirname(script_path)

if len(sys.argv) == 3 and '/' in os.path.abspath(sys.argv[1]):
    setting_file = os.path.abspath(sys.argv[1])
    dir_path = os.path.dirname(setting_file)
    sys.path.insert(0, dir_path)
    from setting import *
else:

    setting_file = ''
    exit('Please using `aft_pipelines_main.py setting.py 1,4,5,6,7,8,9`.The setting.py should in your project path.')
    from setting import *

# relative input csv path
pipelines_path = os.path.join(os.path.dirname(base_outpath),
                              'pipelines_result')

somatic_csv = os.path.join(pipelines_path, 'somatic/')
germline_csv = os.path.join(pipelines_path, 'germline/')
vcf_path = os.path.join(pipelines_path, 'vcf_storge/')
info_summary_path = os.path.join(pipelines_path, 'summary_info/')
# relative output path
final_csv = os.path.join(pipelines_path, 'final_output/')
filtered_csv = os.path.join(pipelines_path, 'filtered_file/')
pcgr_output = os.path.join(pipelines_path, 'pcgr/')
summary_depth = os.path.join(pipelines_path, 'info_summary/')
pack_result = os.path.join(pipelines_path, 'packet_result')


def check_path(path):
    if not os.path.isdir(path):
        os.system('mkdir -p %s ' % path)


for opath in [somatic_csv, germline_csv, vcf_path, final_csv, filtered_csv, pcgr_output, pack_result,
              info_summary_path]:
    check_path(opath)


def run(cmd):
    print(cmd)
    os.system(cmd)


def move_cp(only_info_summary=False):
    if only_info_summary:
        # target cov_info_summary
        cmdline = 'cp %s/%s_result/*/*_cov_summary.info %s' % (
            base_outpath, PROJECT_NAME, info_summary_path)
        run(cmdline)
        exit('Finished download all need file.')
    # target csv file
    cmdline = 'cp %s/%s_result/*_somatic/*anno.csv %s' % (base_outpath, PROJECT_NAME,somatic_csv)
    run(cmdline)
    cmdline = 'cp %s/%s_result/*W/*anno.csv %s' % (base_outpath, PROJECT_NAME, germline_csv)
    run(cmdline)
    # target vcf files
    cmdline = "cp %s/%s_result/*_somatic/*.mt2.vcf %s" % (base_outpath, PROJECT_NAME, vcf_path)
    run(cmdline)
    cmdline = 'cp %s/%s_result/*W/*merged.vcf %s' % (base_outpath, PROJECT_NAME, vcf_path)
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
    cmdline = "cp %s %s" % (fetch_path, final_csv)
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
    file_df.loc[high_cov_bucket, :].to_csv(outpath.replace('.csv', '_f10.csv'), index=False)


def pack_its_up():
    cmdline1 = 'mkdir -p %s/variants_csv; cp %s/*_all_except_AF_PASS_with_info.csv %s/variants_csv/' % (
        pack_result, final_csv, pack_result)
    # cmdline2 = 'mkdir -p %s/html_report; cp %s/*/*.html %s/html_report/' % (pack_result, pcgr_output, pack_result)
    cmdline2 = 'cp %s %s/' % (os.path.join(base_outpath, 'quality_accessment_raw.csv'), pack_result)
    cmdline3 = 'cp %s %s/' % (os.path.join(info_summary_path, 'depth_dis_lines.html'), pack_result)
    cmdline4 = """mkdir -p %s/pre_html_report;; cp %s/*_TS_merged.mt2.vt.vcf %s/pre_html_report/""" % (pack_result,vcf_path,pack_result)
    from simpy_summary import make_summary
    function_csv = make_summary(glob.glob('%s/variants_csv/*.csv' % pack_result))
    function_csv.to_csv('%s/variation_function_stats.csv' % pack_result, index=False)
    run(cmdline1)
    run(cmdline2)
    run(cmdline3)
    run(cmdline4)


def run_cal_cov():
    noexist_cov_bams = []
    for bam in glob.glob("%s/XK_result/*/*sorted.bam" % base_outpath):
        if not os.path.isfile(bam.partition('.')[0] + '_cov.info'):
            noexist_cov_bams.append(bam)
    for bam in glob.glob("%s/XK_result/*/*recal_reads.bam" % base_outpath):
        if not os.path.isfile(bam.partition('.')[0] + '_cov.info'):
            noexist_cov_bams.append(bam)
    if noexist_cov_bams:
        print("There are listed below bam file doesn't cal coverage: %s.\n We need to cal it first. \n" % '\n'.join(
                noexist_cov_bams))
        cmd = """for each in %s; do python %s/pre_pipelines_analysis/cal_Cov_script.py -b $each -B %s -r %s & done""" % (
            ' '.join(noexist_cov_bams),
            os.path.dirname(dir_script.lstrip('/')),
            bed_file_path,
            REF_file_path)
        run(cmd)
    else:
        print('run_cal_cov is fine.')

def run_quality_accessment():
    run_cal_cov()
    if not os.path.isfile(os.path.join(base_outpath, 'quality_accessment_raw.csv')):
        print('Runing quality accessment script.')
        cmd = """/home/liaoth/tools/pysamstats_venv/bin/python2.7 %s/pre_pipelines_analysis/quality_accessment.py %s""" % (
            os.path.dirname(dir_script.lstrip('/')),
            os.path.abspath(
                sys.argv[
                    1]))
        run(cmd)
    else:
        print('Finishing quality accessment script.')

def run_cov_summary():
    run_cal_cov()
    if not glob.glob("%s/XK_result/*/*_summary.info" % base_outpath):
        print('No coverage summary file, we need to cal it first.')
        cmd = """python %s/aft_pipelines_analysis/run_info_summary.py %s""" % (dir_script,os.path.abspath(sys.argv[1]))
        run(cmd)
    else:
        print('run_cov_summary is fine.')

def draw_coverage_depths():
    from draw_quality_line import draw_coverage_depths as dcd
    dcd(info_summary_path, NORMAL_SIG, TUMOR_SIG)
    print('finish drawing.')


if __name__ == '__main__':
    args = sys.argv[-1].split(',')
    print('Using 1,4,5,6,7,8,9 as parameter to init all process. Or it will not do anything.')
    print('Moving all files and processing all files in local path: %s' % (
        pipelines_path))

    server_setting_path = os.path.abspath(sys.argv[1])
    if '1' in args:
        move_cp()
        print('Collecting all needed files into directory.')
        print('Prepareing do cal coverage and quality accessment, it will try to detect old first. If you want to overwrite it. Please delete it by yourself.')
        run_quality_accessment()
    if '2' in args:
        run_cov_summary()
        move_cp(only_info_summary=True)
        print('Collecting all needed files into directory.')
    if '3' in args:
        draw_coverage_depths()
    if '4' in args:
        print('Runing filtered pipelines. Output all samples filtered file into %s' % filtered_csv)
        run_pipelines()
    if '5' in args:
        print("Preparing run 'run_add_per_info_into_csv.py' scripts.")
        run('python %s/run_add_per_info_into_csv.py %s' % (dir_script, server_setting_path))
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
                vt_path=vt_pro, hg19_ref_path=REF_file_path)
            run(cmdline)
        #
        # for vt_vcf in glob.glob('%s/*_TS_merged.mt2.vt.vcf' % vcf_path):
        #     pair_name = os.path.basename(vt_vcf).split('_TS')[0]
        #     cmdline = 'mkdir -p %s/%s_FINAL;' % (pcgr_output, pair_name)
        #     cmdline += "sudo python {pcgr_dir}/pcgr.py --input_vcf {in_vcf} {pcgr_dir} '{pcgr_output}/{pair_name}_FINAL' {pcgr_dir}/data/pcgr_configuration_default.toml {pair_name}_FINAL --force_overwrite"
        #     cmdline = cmdline.format(in_vcf=vt_vcf, pcgr_output=pcgr_output, pair_name=pair_name, pcgr_dir=pcgr_pro)
        #     run(cmdline)
    print("*" * 50)
    print("Finishing all aft pipelines work. Congratulation~~~For Fun!")
    print("*" * 50)
    if '9' in args:
        pack_its_up()
        print("All needed is packed up, you just need to download `pack_result` directory. And run pcgr or otherthing else.")
