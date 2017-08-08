import argparse, sys,os,glob
from luigi.cmdline import luigi_run
from somatic_filter_script import pair_filter,single_filter
import parse_PCR2bed
import result_csv_plus_depth_info

parser_2 = argparse.ArgumentParser()
parser_2.add_argument('-set', '--setting', dest='setting_path', type=str,
                    help="setting path")
parser_2.add_argument('-A', '--analyze', dest='analyze_way', type=str,
                    help="the way you choose to analyze data.[somatic,germline]")
parser_2.add_argument('-cal', '--Cal_cov', dest='cal_cov_info', type=str, help="cal cov info as your setting file.")
parser_2.add_argument('-arg', '--arguments', dest='arguments_space', type=str,
                    help="place you need to type your sample name.")
parser_2.add_argument('-parse', '--parsePCR', dest='PCR_PATH', type=str, help="Enter blastn result path.")
parser_2.add_argument('-fix', dest='fix_format', type=str, help="formatted for csv result file.")
parser_2.add_argument('-SBF', '--Somatic_combined_file', dest='SBF', type=bool,
                    help="Turn on/off the single analysis after SomaticPipelines.")
parser_2.add_argument('-SSF', '--Somatic_single_file', dest='SSF', type=bool,
                    help="Turn on/off the single analysis after SomaticPipelines.")

args = parser_2.parse_args()


setting_path = args.setting_path
SBF = args.SBF
SSF = args.SSF
analyze_way = args.analyze_way
cal_cov_info = args.cal_cov_info
arguments_space = args.arguments_space
PCR_PATH = args.PCR_PATH
fix_format = args.fix_format

if setting_path:
    import_path = setting_path.rpartition('/')[0]
    #import_name = setting_path.rpartition('/')[2][:-3]
    sys.path.insert(0, import_path)
    from setting import *
    from parse_file_name import pfn
else:
    from setting import *
    from parse_file_name import pfn

def formatter_output(args):
    if ',' in args:
        arg_List = args.split(',')
        parsed_name = []
        for val in arg_List:
            parsed_name.append(pfn(PE1_fmt.format(input=val), 'all'))
    else:
        val = args
        parsed_name = pfn(PE1_fmt.format(input=val),'all')

    if not self_adjust_fn:
        fq_file = PE1_fmt.format(input=val)
    else:
        input_list = glob.glob(base_inpath + '/*' + val + '*')
        if filter_str:
            input_list = [_i.replace(fq_suffix, '') for _i in input_list if filter_str not in _i]
        fq_file = input_list[0]

    output = """Current Variants: Please make sure your variants is right.\n\n
    Input path: {b_i}
    output path: {b_o}
    Sig represent NORMAL: {sig_n}
    Sig represent TUMOR: {sig_T}
    Pair file format: {pe_fmt}
    input_fastq_file example:{fq}\n

    One of args filename parsed: {parsed_result}""".format(b_i=base_inpath, b_o=base_outpath, sig_n=NORMAL_SIG,
                                                           sig_T=TUMOR_SIG, pe_fmt=PE1_fmt, fq=fq_file,
                                                           parsed_result=str(parsed_name))
    return output

if __name__ == '__main__':

    GO_ON = True
    if SBF:
        inpath = raw_input('Please Enter the base input path \n')
        outpath = raw_input('Please Enter the base output path \n')
        sample_name = raw_input('Please Enter the sample name without normal/tumor sig. Like XK-2.\n')
        _input1 = raw_input('Please Enter the filename which is normal only analysis from Mutect2 (.csv)\n')
        _input2 = raw_input('Please Enter the filename which is tumor only analysis from Mutect2 (.csv)\n')
        input1 = inpath + '/' + _input1
        input2 = inpath + '/' + _input2
        output1 = outpath + '/' + sample_name + '_mt2_low_F.csv'
        output2 = outpath + '/' + sample_name + '_mt2_NotInNormal.csv'
        pair_filter(input1,input2,output1,output2)
    if SSF:
        inpath = raw_input('Please Enter the base input path \n')
        outpath = raw_input('Please Enter the base output path \n')
        sample_name = raw_input('Please Enter the sample name without normal/tumor sig. Like XK-2.\n')
        _input1 = raw_input('Please Enter the filename which is pair analysis from Mutect2 (.csv)\n')
        input1 = inpath + '/' + _input1
        output1 = outpath + '/' + sample_name + '_mt2_filted.csv'
        single_filter(input1,output1,disease= False)
    if PCR_PATH and GO_ON:
        op_path = PCR_PATH.rpartition('.')[0]+'.bed'
        if os.path.isfile(op_path):
            print 'Bed file already exist.'
        else:
            print 'Path of pcr blastn file: '+ PCR_PATH +'\n'
            print 'Total gene: ' + str(total_gen) +'\n'
            print 'Output path of bed file: ' + op_path +'\n'
            parse_PCR2bed.parse_PCR2BED(PCR_PATH,op_path,total_gen)
            print 'PCR2bed completed.'

    if arguments_space and analyze_way and GO_ON:
        reload(sys)
        parse_args = arguments_space.split(',')
        if '.fastq.gz' in parse_args[0]:
            parse_args = [i.rpartition('.fastq.gz')[0] for i in parse_args]
        parse_args = ','.join(parse_args)
        sys.path.append(os.path.abspath(sys.argv[0]))
        if analyze_way == 'somatic':
            RUN_COMMAND = raw_input("""%s\n\n If you sure, please input Y/y.""" % formatter_output(arguments_space))
            if RUN_COMMAND.upper() == 'Y':
                cmd_str = "--module SomaticPipelines workflow --x {args} --local-scheduler --workers {worker}".format(args=parse_args, worker=str(worker))
                luigi_run(cmd_str.split(' '))
            else:
                print 'Exit now.'
                GO_ON=False
        elif analyze_way == 'germline':
            RUN_COMMAND = raw_input("""%s\n\n If you sure, please input Y/y.""" % formatter_output(arguments_space))
            if RUN_COMMAND.upper() == 'Y':
                cmd_str ='--module GermlinePipelines workflow --x {args} --local-scheduler --workers {worker}'.format(args=parse_args, worker=str(worker))
                luigi_run(cmd_str.split(' '))
            else:
                print 'Exit now.'
                GO_ON=False

        elif analyze_way == 'somatic_gemini':
            RUN_COMMAND = raw_input("""%s\n\n If you sure, please input Y/y.""" % formatter_output(arguments_space))
            if RUN_COMMAND.upper() == 'Y':
                cmd_str ="--module SomaticPipelines_to_gemini workflow --x {args} --local-scheduler --workers {worker}".format(args=parse_args, worker=str(worker))
                luigi_run(cmd_str.split(' '))
            else:
                print 'Exit now.'
                GO_ON=False

        elif analyze_way == 'analyze_gemini':
            RUN_COMMAND = raw_input("""%s\n\n If you sure, please input Y/y.""" % formatter_output(arguments_space))
            if RUN_COMMAND.upper() == 'Y':
                cmd_str = "--module GermlinePipelines_to_gemini workflow --x {args} --local-scheduler --workers {worker}".format(args=parse_args, worker=str(worker))
                luigi_run(cmd_str.split(' '))
            else:
                print 'Exit now.'
                GO_ON=False
        elif analyze_way == 'test_multi':
            RUN_COMMAND = raw_input("""%s\n\n If you sure, please input Y/y.""" % formatter_output(arguments_space))
            if RUN_COMMAND.upper() == 'Y':
                cmd_str = "--module test_multi_SomaticPipelines workflow --x {args} --local-scheduler --workers {worker}".format(args=parse_args, worker=str(worker))
                luigi_run(cmd_str.split(' '))
            else:
                print 'Exit now.'
                GO_ON=False

        else:
            print 'Please use correct args in -A, either somatic/germline'
    else:
        print 'Fatal error occur, please enter both -A and -arg.'

    if cal_cov_info and GO_ON:
        print 'Using bed file : '+bed_file_path
        print 'Using bam file : '+ output_fmt.format(path=cal_cov_info, SN=cal_sample_name)+'.recal_reads.bam'
        print 'Output file will be output at : ' + output_fmt.format(path=cal_cov_info, SN=cal_sample_name)+'_cov.info'
        cal_sample_name(cal_cov_info,cal_sample_name,bed_file_path)
        print 'cal_cov_info Completed.'
    else:
        print "Maybe the file you want to cal doesn't exist."

    if fix_format and GO_ON:
        try:
            depth_info,csv_result = fix_format.split(',')
        except:
            print "You need to type args like this, 'depth_info_path','csv_result_path'"
        if not os.path.isfile(depth_info):
            print "depth_info file doesn't exist"
        if not os.path.isfile(csv_result):
            print "csv_result file doesn't exist"

        print 'Fixing file : ' + csv_result + '  Using depth file : ' + depth_info
        fixed_file = result_csv_plus_depth_info.plus_depth_info(depth_info,csv_result)
        with open(csv_result.rpartition('.csv')[0]+'formatted.csv','w') as f1:
            fixed_file.to_csv(f1,index=False)

        # print """e.g
        # 1. When you want to do somatic/germline analysis.Like this:
        #     python main.py -s /home/liaoth/project_formal/170123_XK/setting.py -A somatic -arg XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18
        # 2. When you analyze somatic data, you need to put it into a batch process program.
        #     python main.py -
        # """

