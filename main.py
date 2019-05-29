import argparse
import glob
import os
import re
import sys

import parse_PCR2bed
from parse_file_name import fileparser
from aft_pipelines_analysis import result_csv_plus_depth_info
from luigi.cmdline import luigi_run
from somatic_filter_script import pair_filter,single_filter


argparser = argparse.ArgumentParser()

argparser.add_argument('-i', '--input', dest='input_file', type=str,
                    help="data tab for process")
argparser.add_argument('-set', '--setting', dest='setting_path', type=str,
                    help="setting path")
argparser.add_argument('-A', '--analyze', dest='analyze_way', type=str,
                    help="the way you choose to analyze data.[somatic,germline,germline_gemini,somatic_gemini]")
argparser.add_argument('-cal', '--Cal_cov', dest='cal_cov_info', action="store_true",
                      help="Cal cov info for quality assessment file.")
argparser.add_argument('-arg', '--arguments', dest='arguments_space', type=str,
                    help="place you need to type your sample name.(if you want to use lots of in your dir,please type in 'auto_detect'.)")
argparser.add_argument('-parse', '--parsePCR', dest='PCR_PATH', type=str, help="Enter blastn result path.")
argparser.add_argument('-fix', dest='fix_format', type=str, help="formatted for csv result file.")
argparser.add_argument('-SBF', '--Somatic_combined_file', dest='SBF', type=bool,
                    help="Turn on/off the single analysis after SomaticPipelines.")
argparser.add_argument('-SSF', '--Somatic_single_file', dest='SSF', type=bool,
                    help="Turn on/off the single analysis after SomaticPipelines.")
argparser.add_argument('-debug', '--debug', dest='debug_', action="store_true", help="Entering debug mode.")
args = argparser.parse_args()

data_csv = args.input_file
setting_path = args.setting_path
SBF = args.SBF
SSF = args.SSF
analyze_way = args.analyze_way
cal_cov_info = args.cal_cov_info
arguments_space = args.arguments_space
PCR_PATH = args.PCR_PATH
fix_format = args.fix_format
debug_ = args.debug_
if setting_path:
    import_path = os.path.dirname(setting_path)
    import_name = os.path.basename(setting_path).replace('.py', '')
    sys.path.insert(0, import_path)
    __import__(import_name, globals=globals())
else:
    from setting import *


if __name__ == '__main__':
    # if SBF:
    #     inpath = input('Please Enter the base input path \n')
    #     outpath = input('Please Enter the base output path \n')
    #     sample_name = input('Please Enter the sample name without normal/tumor sig. Like XK-2.\n')
    #     _input1 = input('Please Enter the filename which is normal only analysis from Mutect2 (.csv)\n')
    #     _input2 = input('Please Enter the filename which is tumor only analysis from Mutect2 (.csv)\n')
    #     input1 = inpath + '/' + _input1
    #     input2 = inpath + '/' + _input2
    #     output1 = outpath + '/' + sample_name + '_mt2_low_F.csv'
    #     output2 = outpath + '/' + sample_name + '_mt2_NotInNormal.csv'
    #     pair_filter(input1,input2,output1,output2)
    # if SSF:
    #     inpath = input('Please Enter the base input path \n')
    #     outpath = input('Please Enter the base output path \n')
    #     sample_name = input('Please Enter the sample name without normal/tumor sig. Like XK-2.\n')
    #     _input1 = input('Please Enter the filename which is pair analysis from Mutect2 (.csv)\n')
    #     input1 = inpath + '/' + _input1
    #     output1 = outpath + '/' + sample_name + '_mt2_filted.csv'
    #     single_filter(input1,output1,disease= False)

    if PCR_PATH:
        # todo: table fit it or not?
        op_path = PCR_PATH.rpartition('.')[0]+'.bed'
        if os.path.isfile(op_path):
            print('Bed file already exist.')
        else:
            print('Path of pcr blastn file: '+ PCR_PATH +'\n')
            print('Total gene: ' + str(total_gen) +'\n')
            print('Output path of bed file: ' + op_path +'\n')
            parse_PCR2bed.parse_PCR2BED(PCR_PATH,op_path,total_gen)
            print('PCR2bed completed.')

    if arguments_space and analyze_way:

        sys.path.append(os.path.abspath(sys.argv[0]))
        sys.path.append(os.path.dirname(__file__)+'/luigi_pipelines/')

        if bip:
            scheduler_str = '--scheduler-host {bip}'.format(bip=bip)
        else:
            scheduler_str = '--local-scheduler'

        if analyze_way == 'somatic':

            cmd_str = "--module SomaticPipelines workflow --x {args} {sch} --workers {worker}".format(args=parse_args,
                                                                                                      worker=str(worker),
                                                                                                      sch=scheduler_str)
            luigi_run(cmd_str.split(' '))

        elif analyze_way == 'germline':


            cmd_str ='--module GermlinePipelines workflow --x {args} {sch} --workers {worker}'.format(args=parse_args, worker=str(worker),sch=scheduler_str)
            luigi_run(cmd_str.split(' '))

        elif analyze_way == 'somatic_gemini':
            cmd_str ="--module SomaticPipelines_to_gemini workflow --x {args} {sch} --workers {worker}".format(args=parse_args, worker=str(worker),sch=scheduler_str)
            luigi_run(cmd_str.split(' '))

        elif analyze_way == 'germline_gemini':
            cmd_str = "--module GermlinePipelines_to_gemini workflow --x {args} {sch} --workers {worker}".format(args=parse_args, worker=str(worker),sch=scheduler_str)
            luigi_run(cmd_str.split(' '))

        elif analyze_way == 'somatic_gatk4':
            cmd_str = "--module SomaticPipelines_gatk4 workflow --x {args} {sch} --workers {worker}".format(
                args=parse_args, worker=str(worker), sch=scheduler_str)
            luigi_run(cmd_str.split(' '))

        elif analyze_way == 'germline_gatk4':
            cmd_str = "--module GermlinePipelines_gatk4 workflow --x {args} {sch} --workers {worker}".format(
                args=parse_args, worker=str(worker), sch=scheduler_str)
            luigi_run(cmd_str.split(' '))
        else:
            exit('Please use correct args in -A, either somatic/germline')
    else:
        exit('Fatal error occur, please enter both -A and -arg.')

    if cal_cov_info:
        cmdline = 'python %s %s' % (
        '/home/liaoth/project/Whole_pipelines/pre_pipelines_analysis/quality_accessment.py', setting_path)
        os.system(cmdline)
    #
    # if fix_format:
    #     try:
    #         depth_info,csv_result = fix_format.split(',')
    #     except:
    #         print("You need to type args like this, 'depth_info_path','csv_result_path'")
    #     if not os.path.isfile(depth_info):
    #         print("depth_info file doesn't exist")
    #     if not os.path.isfile(csv_result):
    #         print("csv_result file doesn't exist")
    #
    #     print('Fixing file : ' + csv_result + '  Using depth file : ' + depth_info)
    #     fixed_file = result_csv_plus_depth_info.plus_depth_info(depth_info, csv_result)
    #     with open(csv_result.rpartition('.csv')[0]+'formatted.csv','w') as f1:
    #         fixed_file.to_csv(f1,index=False)

# print("""e.g)
# 1. When you want to do somatic/germline analysis.Like this:
#     python main.py -s /home/liaoth/project_formal/170123_XK/setting.py -A somatic -arg XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18
# 2. When you analyze somatic data, you need to put it into a batch process program.
#     python main.py -
# """
