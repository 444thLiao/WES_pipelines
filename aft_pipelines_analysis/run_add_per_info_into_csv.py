
import glob
import os,sys
import multiprocessing

script_path = __file__
############################################################
python_exe = '/home/liaoth/tools/pysamstats_venv/bin/python'

cmdline = '{py_exe} %s/add_per_info_into_csv.py -i {input} -o {output} -tb {tb} -nb {nb} -b /home/liaoth/data_bank/XK_WES/Sureselect_V6_COSMIC_formal.bed' % os.path.dirname(script_path)


def run(cmd):
    print(cmd)
    os.system(cmd)

def run_batch(input_path,tb_bam_path,nb_bam_path,NORMAL_SIG,TUMOR_SIG):
    cmdlines = []
    for path in glob.glob(input_path):
        each_pair = os.path.basename(path).split('_all_except')[0]
        if each_pair.count('-') == 2:
            normal_name = each_pair.rpartition('-')[0] + NORMAL_SIG
            tumor_name = each_pair.rpartition('-')[0] + TUMOR_SIG + '-' + each_pair.rpartition('-')[2]
        else:
            normal_name = each_pair + NORMAL_SIG
            tumor_name = each_pair + TUMOR_SIG
        output = output_path.format(prefix=path.split('_all_except')[0])
        tb = tb_bam_path.format(tb=tumor_name)
        nb = nb_bam_path.format(nb=normal_name)
        cmdlines.append(cmdline.format(py_exe=python_exe,input=path,output=output,tb=tb,nb=nb))
    return cmdlines
if __name__ == '__main__':
    from ..setting import *
    if len(sys.argv) ==2 and '/' in sys.argv[-1]:
        setting_file = os.path.abspath(sys.argv[-1])
        dir_path = os.path.dirname(setting_file)
        sys.path.insert(0,dir_path)
        exec('from setting import *')
    else:
        print('Please using `run_add_per_info_into_csv.py setting.py`.The setting.py should in your project path.')
        exit()

    num_processes = 4
    input_path = '%s/temp_/*_all_except_AF_depth_PASS.csv' % base_inpath
    output_path = '{prefix}_all_except_AF_PASS_with_info.csv'
    tb_bam_path = "%s/output/XK_result/{tb}/{tb}.recal_reads.bam" % base_inpath
    nb_bam_path = "%s/output/XK_result/{nb}/{nb}.recal_reads.bam" % base_inpath

    cmdlines = run_batch(input_path,tb_bam_path,nb_bam_path,NORMAL_SIG,TUMOR_SIG)
    makesure = input("If your `num of processes >4`, Please be careful of memory. It may stalled whole server.\nUsing %s processes, prepare process listing files: \n . %sIf you make sure, please type y/Y." % (num_processes,'\n'.join(glob.glob(input_path))))

    if makesure.strip().upper() == 'Y':
        pool = multiprocessing.Pool(num_processes)
        pool.map(run, cmdlines)
    else:
        print('Exiting ......')
        exit()
