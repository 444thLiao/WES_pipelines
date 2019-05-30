import glob
import os
import threading
from tqdm import tqdm

# obsoleted file

def format_cov_depth_heatmap(info_file, genes, block=15):
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
        for _b in range(0, block):
            block_start = total_len / block * _b
            block_next = total_len / block * (_b + 1)
            if _b != block - 1:
                block_df = piece_df.iloc[block_start:block_next]
            else:
                block_df = piece_df.iloc[block_start:]
            if len(block_df) == 0:
                init_blocks_ave_depth.append(0)
                continue
            block_total_depth = block_df.coverage.sum()
            block_ave_depth = block_total_depth / float(len(block_df))
            init_blocks_ave_depth.append(block_ave_depth)
        result_df.loc[_gene, :] = init_blocks_ave_depth
    result_df.to_csv(info_file.replace('_cov.info', '_cov_gene.heatmap'))
    return result_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-input_pattern', dest='pattern1', required=True, type=str)
    parser.add_argument('-name_patern', dest='pattern2', required=True, type=str)
    parser.add_argument('-cache', dest='cache', required=True, type=str)
    parser.add_argument('-output', dest='output', required=True, type=str)

    args = parser.parse_args()
    pattern = args.pattern
    pattern2 = args.pattern2
    output = args.output
    cache = args.cache

    for path in glob.glob('/home/liaoth/project/180104_XK/output/XK_result/*/*_cov.info'):
        if 'sorted' not in path and not os.path.isfile(path.replace('_cov.info', '_cov_summary.info')):
            t = threading.Thread(target=cov_depth_to_file, kwargs={'cov_info': path})
            print('processing......', path)
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
            print('processing......', path)
            t.setDaemon(True)
            t.start()
