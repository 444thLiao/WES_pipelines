import re, glob, pandas, os
from pandas import DataFrame as df
import time, tqdm
import threading
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

t1 = time.clock()
tmp = cov_depth('~/project/XK_WES/180309_all/output/XK_result/XK-25T/XK-25T_cov.info')
print(time.clock()-t1)

def cov_depth_to_file(cov_info, output_each_pos=False):
    """
    cal the info from cal_Cov.py and produce a plot
    :param cov_info:
    :return:
    """
    tmp = open(cov_info).xreadlines()
    total_lines = int(os.popen('wc -l %s' % cov_info).read().split(' ')[0])

    bases = ['A', 'T', 'C', 'G']
    each_pos = []
    pos_list = []
    gene_list = []
    chr_list = []
    count = 0
    pbar = ProgressBar()
    pbar = pbar.start()
    for line in tmp:
        if count == 0:
            idx = [line.split('\t').index(_) for _ in bases]
        else:
            each_pos.append(sum([int(line.split('\t')[_]) for _ in idx]))
            pos_list.append(line.split('\t')[2])
            gene_list.append(line.split('\t')[0].replace('ref|', ''))
            chr_list.append(line.split('\t')[1])
        count += 1
        pbar.update(int((float(count) / (total_lines)) * 100))
    if output_each_pos:
        parsed_df = df.from_dict(dict(pos=pos_list, coverage=each_pos, gene_name=gene_list, chr_name=chr_list))
        return parsed_df
    pbar.finish()
    result_depths = []
    result_coverages = []
    depths_name = []
    for _i in tqdm.tqdm(range(1, 101)):
        _coverage = len([_ for _ in each_pos if _ >= _i * 10])
        _depth = _i * 10
        depths_name.append('%sX' % str(_depth))
        result_coverages.append(_coverage)
        result_depths.append(_depth)
    df_result = df(index=depths_name, columns=['coverage', 'depths'])
    df_result.loc[:, 'coverage'] = result_coverages
    df_result.loc[:, 'depths'] = result_depths
    df_result.to_csv(cov_info.replace('_cov.info', '_cov_summary.info'), sep='\t')

def format_cov_depth_heatmap(info_file,genes,block=15):
    if type(info_file) == str:
        info_df = cov_depth_to_file(info_file, output_each_pos=True)
    else:
        info_df = info_file
    result_df = df(index=genes, columns=range(block))
    for _gene in tqdm.tqdm(genes):
        piece_df = info_df[info_df.loc[:, 'gene_name'] == _gene]
        total_len = len(piece_df)
        piece_df = piece_df.sort_values('pos')

        init_blocks_ave_depth = []
        for _b in xrange(0,block):
            block_start = total_len/block * _b
            block_next = total_len / block * (_b+1)
            if _b != block-1:
                block_df = piece_df.iloc[block_start:block_next]
            else:
                block_df = piece_df.iloc[block_start:]
            if len(block_df) == 0:
                init_blocks_ave_depth.append(0)
                continue
            block_total_depth = block_df.coverage.sum()
            block_ave_depth = block_total_depth/float(len(block_df))
            init_blocks_ave_depth.append(block_ave_depth)
        result_df.loc[_gene, :] = init_blocks_ave_depth
    result_df.to_csv(info_file.replace('_cov.info', '_cov_gene.heatmap'))
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

    for path in glob.glob('/home/liaoth/project/180104_XK/output/XK_result/*/*_cov.info'):
        if 'sorted' not in path and not os.path.isfile(path.replace('_cov.info', '_cov_summary.info')):
            t = threading.Thread(target=cov_depth_to_file, kwargs={'cov_info': path})
            print 'processing......', path
            t.setDaemon(True)
            t.start()
    sum_sig_genes = ['AKT1', 'AKT2', 'AKT3', 'ARID1A', 'ARID1B', 'ARID2', 'ASCL4', 'ATM', 'BRAF', 'CDKN2A', 'COBL',
                     'CREBBP', 'CTNNB1', 'CUL3', 'EGFR', 'EP300', 'EPHA7', 'ERBB2', 'ERBB3', 'FGFR1', 'FGFR2', 'FGFR3',
                     'FOXP2', 'HRAS', 'KEAP1', 'KMT2D', 'KRAS', 'MAP2K1', 'MET', 'MGA', 'KMT2A', 'NF1', 'NFE2L2',
                     'NOTCH1',
                     'NOTCH2', 'NRAS', 'PIK3CA', 'PTEN', 'RASA1', 'RB1', 'RBM10', 'RIT1', 'SETD2', 'SLIT2', 'SMAD4',
                     'SMARCA4', 'SOX2-OT', 'STK11', 'TP53', 'TP63', 'TSC1', 'TSC2', 'U2AF1', 'ALK', 'DDR2']

    # format_cov_depth_heatmap('/home/liaoth/project/180104_XK/output/XK_result/XK-32W/XK-32W_cov.info',sum_sig_genes)
    for path in glob.glob('/home/liaoth/project/180104_XK/output/XK_result/*/*_cov.info'):
        if 'sorted' not in path and not os.path.isfile(path.replace('_cov.info', '_cov_gene.heatmap')):
            t = threading.Thread(target=format_cov_depth_heatmap, kwargs={'info_file': path, 'genes': sum_sig_genes})
            t.setName('deal with %s' % os.path.basename(path).split('_cov')[0])
            print 'processing......', path
            t.setDaemon(True)
            t.start()
