import re, glob, pandas, os
from pandas import DataFrame as df
import time, tqdm
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
    for _i in xrange(1,101):
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


pool = multiprocessing.Pool(4)
pool.map(cov_depth,[_ for _ in glob.glob('/home/liaoth/project/XK_WES/180309_all/output/XK_result/*/*cov.info') if 'sorted' not in _])
