import re,glob,pandas
from pandas import DataFrame as df

import time
def cov_depth(cov_info):
    """
    cal the info from cal_Cov.py and produce a plot
    :param cov_info:
    :return:
    """
    raw_info = pandas.read_csv(cov_info,sep='\t',index_col=False,engine='python')
    raw_info.base = list(raw_info.loc[:, ['A', 'T', 'C', 'G']].sum(1))

    #total = len(raw_info)
    #coverage_info = list(raw_info.loc[:, ['A', 'T', 'C', 'G']].sum(1))

    result_depths = []
    result_coverages = []
    depths_name = []
    for _i in xrange(1,101):
        _coverage = len(raw_info[raw_info.base >= _i*10])
        _depth = _i*10
        depths_name.append('%sX' % str(_depth))
        result_coverages.append(_coverage)
        result_depths.append(_depth)
    return result_depths,result_coverages,depths_name


def format_cov_depth_heatmap(info_file,genes,block=15):
    info_df = df.from_csv(info_file,sep='\t',index_col=False)

    result_df = df(index=genes,columns=range(block))
    for _gene in genes:
        piece_df = info_df[info_df.Gene == _gene]
        total_len = len(piece_df)
        piece_df = piece_df.sort_values('Posistion')

        init_blocks_ave_depth = []
        for _b in xrange(0,block):
            block_start = total_len/block * _b
            block_next = total_len / block * (_b+1)
            if _b != block-1:
                block_df = piece_df.iloc[block_start:block_next]
            else:
                block_df = piece_df.iloc[block_start:]
            block_total_depth = sum(block_df.sum(0).loc[['A', 'T', 'C', 'G']])
            block_ave_depth = block_total_depth/float(len(block_df))
            init_blocks_ave_depth.append(block_ave_depth)
        result_df.loc[_gene] = init_blocks_ave_depth
    return result_df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_pattern',dest='pattern1',required=True,type=str)
    parser.add_argument('-name_patern',dest='pattern2',required=True,type=str)
    parser.add_argument('-cache', dest='cache', required=True, type=str)
    parser.add_argument('-output',dest='output',required=True,type=str)

    args = parser.parse_args()
    pattern = args.pattern
    pattern2=args.pattern2
    output=args.output
    cache=args.cache
    df_bucket = []
    #for path in glob.glob('/home/liaoth/project/YZ_XR_WES/output/YZ_XR_WES_result/HCC591*/HCC591*_cov.info'):
    for path in glob.glob('/media/vd1/past_project/170123_XK/output/XK_result/*/*_cov.info'):
    #for path in glob.glob(pattern):
        print 'processing......',path
        t1=time.time()
        depths,coverage,depths_name = cov_depth(path)
        _cache = df(columns=['dep','cov','sample'])
        _cache.dep = depths
        _cache.iloc[:,1] = coverage
        # input_pair_tuple = [('XR7B9017', 'XR7F9013'),
        #                     ('XR7B9014', 'XR7B9013'),
        #                     ('XR7B9012', 'XR7B9011'),
        #                     ('XR7B9010', 'XR7B9009')]

        #sample_name = re.findall(pattern2,path)[0]
        sample_name = re.findall('(XK-[28][TW])_S[0-9]{2}_', path)[0]
        _cache.iloc[:, 2] = sample_name

        df_bucket.append(_cache)
        print 'done...one.....used',time.time()-t1

    a = pandas.concat(df_bucket)
    with open(cache,'w') as f1:
        a.to_csv(f1)

    #
    import seaborn as sns
    ax = sns.stripplot(x='dep',y='cov',hue='sample',data=a)

    ######   WES all base count = 60456963

    depths_name = []
    for _i in xrange(1,101):
        _depth = _i * 10
        depths_name.append('%sX' % str(_depth))
    depth_name2 = []
    for _idx,_v in enumerate(depths_name):
        if int(_v.replace('X','')) % 50 ==0:
            depth_name2.append(_v)
        else:
            depth_name2.append('')
    ax.set_xticklabels(depth_name2)
    ax.figure.savefig(output)