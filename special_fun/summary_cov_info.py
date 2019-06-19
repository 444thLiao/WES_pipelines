#################################################################################
#### Convert _cov.info into a summary file for plot.
#### validate at 190619 by TianHua Liao
#################################################################################

import numpy as np
import pandas as pd
from tqdm import tqdm
import os

def summarize_covinfo(cov_info, output_f=None):
    """
    cal the info from cal_Cov.py and produce a plot
    :param cov_info:
    :return:
    """
    with open(cov_info, 'r') as f1:
        header_str = next(f1)
        # extract the header and close it...
    header = header_str.strip('\n').split('\t')
    depths_threshold = np.array([10 * _ for _ in range(1, 101)])
    depths_count = np.zeros(len(depths_threshold))
    with open(cov_info, 'r') as f1:
        for row in tqdm(f1):
            if row == header_str:
                continue
            _cache = row.strip('\n').split('\t')
            _cache = dict(zip(header, _cache))
            depth = float(_cache["base"])
            depths_count[depth >= depths_threshold] += 1

    df_result = pd.DataFrame(data=np.array([depths_threshold, depths_count]).T,
                             index=["%sX" % str(int(_)) for _ in depths_threshold],
                             columns=['depths',
                                      'number of pos larger than this depth'])
    if output_f is not None:
        df_result.to_csv(output_f,index=1)
    return df_result

if __name__ == '__main__':
    import sys
    files = sys.argv[1:]
    for f in files:
        abspath_f = os.path.abspath(f)
        summarize_covinfo(abspath_f,
                          abspath_f.replace('cov.info',
                                            'cov_summary.info'))