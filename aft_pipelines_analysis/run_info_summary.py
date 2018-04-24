import re, glob, pandas, os,argparse
from pandas import DataFrame as df
import tqdm,sys
import multiprocessing
def run(cmd):
    print(cmd)
    os.system(cmd)

def cov_depth(cov_info):
    """
    cal the info from cal_Cov.py and produce a plot
    :param cov_info:
    :return:
    """
    raw_info = pandas.read_csv(cov_info,sep='\t',index_col=False,engine='python')
    raw_info.base = list(raw_info.loc[:, ['A', 'T', 'C', 'G']].sum(1))

    result_depths = []
    result_coverages = []
    depths_name = []
    for _i in tqdm.tqdm(range(1,101),total=100):
        _coverage = len(raw_info[raw_info.base >= _i*10])
        _depth = _i*10
        depths_name.append('%sX' % str(_depth))
        result_coverages.append(_coverage)
        result_depths.append(_depth)
    df_result = df(index=depths_name, columns=['coverage', 'depths'])
    df_result.loc[:, 'coverage'] = result_coverages
    df_result.loc[:, 'depths'] = result_depths
    df_result.to_csv(cov_info.replace('_cov.info', '_cov_summary.info'), sep='\t')
    return df_result
    # return result_depths,result_coverages,depths_name

if __name__ == '__main__':
    if len(sys.argv) ==2 and '/' in sys.argv[-1]:
        setting_file = os.path.abspath(sys.argv[-1])
        dir_path = os.path.dirname(setting_file)
        sys.path.insert(0,dir_path)
        from setting import *
    else:
        exit('Please using `run_info_summary.py setting.py`.The setting.py should in your project path.')

    num_processes = 4
    parsing_path = '%s/XK_result/*/*cov.info' % base_outpath

    makesure = str(input("If your `num of processes >4`, Please be careful of memory. It may stalled whole server.\nUsing %s processes, prepare process listing files: \n . %s\n\nIf you make sure, please type y/Y." %
                         (num_processes,'\n'.join(glob.glob(parsing_path)))
                         )
                   )

    if makesure.strip().upper() == 'Y':
        pool = multiprocessing.Pool(num_processes)
        pool.map(cov_depth,[_ for _ in glob.glob(parsing_path) if 'sorted' not in _])
    else:
        print('Exiting ......')
        exit()
