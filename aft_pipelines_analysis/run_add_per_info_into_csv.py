
import glob
import os
import multiprocessing

rep_name = '180321'
############################################################
python_exe = '/home/liaoth/tools/pysamstats_venv/bin/python'

cmdline = '{py_exe} ~/tools/Whole_pipelines/aft_pipelines_analysis/add_per_info_into_csv.py -i {input} -o {output} -tb {tb} -nb {nb} -b ~/data_bank/XK_WES/Sureselect_V6_COSMIC_formal.bed'

input_path = '/home/liaoth/project/XK_WES/%s/temp_/*_all_except_AF_depth_PASS.csv' % rep_name
output_path = '{prefix}_all_except_AF_PASS_with_info.csv'

tb_bam_path = "~/project/XK_WES/%s/output/XK_result/{tb}/{tb}.recal_reads.bam" % rep_name
nb_bam_path = "~/project/XK_WES/%s/output/XK_result/{nb}/{nb}.recal_reads.bam" % rep_name
cmdlines = []
def run(cmd):
    print(cmd)
    os.system(cmd)

for path in glob.glob(input_path):
    each_pair = os.path.basename(path).split('_all_except')[0]
    if each_pair.count('-') == 2:
        normal_name = each_pair.rpartition('-')[0] + 'W'
        tumor_name = each_pair.rpartition('-')[0] + 'T' + '-' + each_pair.rpartition('-')[2]
    else:
        normal_name = each_pair + 'W'
        tumor_name = each_pair + 'T'

    output = output_path.format(prefix=path.split('_all_except')[0])
    tb = tb_bam_path.format(tb=tumor_name)
    nb = nb_bam_path.format(nb=normal_name)

    cmdlines.append(cmdline.format(py_exe=python_exe,input=path,output=output,tb=tb,nb=nb))

pool = multiprocessing.Pool(15)
pool.map(run, cmdlines)
