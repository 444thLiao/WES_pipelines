from pandas import DataFrame as df
import matplotlib.pyplot as plt
import seaborn as sns
import glob,os
import re
sum_sig_genes = ['AKT1', 'AKT2', 'AKT3', 'ARID1A', 'ARID1B', 'ARID2', 'ASCL4', 'ATM', 'BRAF', 'CDKN2A', 'COBL',
                 'CREBBP', 'CTNNB1', 'CUL3', 'EGFR', 'EP300', 'EPHA7', 'ERBB2', 'ERBB3', 'FGFR1', 'FGFR2', 'FGFR3',
                 'FOXP2', 'HRAS', 'KEAP1', 'KMT2D', 'KRAS', 'MAP2K1', 'MET', 'MGA', 'KMT2A', 'NF1', 'NFE2L2', 'NOTCH1',
                 'NOTCH2', 'NRAS', 'PIK3CA', 'PTEN', 'RASA1', 'RB1', 'RBM10', 'RIT1', 'SETD2', 'SLIT2', 'SMAD4',
                 'SMARCA4', 'SOX2-OT', 'STK11', 'TP53', 'TP63', 'TSC1', 'TSC2', 'U2AF1','ALK','DDR2',
                 'CD274','PDCD1LG2','LGALS9']
def extracted_info(info_file,genes,output_dir,suffix=''):
    info_df = df.from_csv(info_file,sep='\t',index_col=False)
    if len(info_df.columns) == 1:
        info_df = df.from_csv(info_file)
    try:
        new_df = info_df[info_df.Gene.isin(genes)]
        if len(new_df) == 0:
            raise IOError
    except:
        new_df = info_df[info_df.Gene.isin(['ref|'+_i for _i in genes])]
    new_df.to_csv(os.path.join(output_dir,os.path.basename(info_file)+suffix))



def format_cov_depth_heatmap(info_file,genes,normalize_coeff,block=15):
    info_df = df.from_csv(info_file,sep='\t',index_col=False)
    if len(info_df.columns) == 1:
        info_df = df.from_csv(info_file)
    result_df = df(index=set(genes),columns=range(block))
    for _gene in set(genes):
        piece_df = info_df[info_df.Gene == 'ref|'+_gene]
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
            block_ave_depth = block_total_depth/float(len(block_df))/normalize_coeff
            init_blocks_ave_depth.append(block_ave_depth)
        try:
            result_df.loc[_gene] = init_blocks_ave_depth
        except:
            import pdb;pdb.set_trace()
    return result_df


import pandas

normalize_df = pandas.read_excel('/home/liaoth/project/170602_XK/WES两次活检总体评估表_改正2_8顺序.xlsx')
count = 1
if __name__ == '__main__':
    for path in [
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-2W_cov_extracted.info',
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-2T_S20_cov_extracted.info',
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-2T-2_cov_extracted.info',
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-8W_cov_extracted.info',
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-8T_S21_cov_extracted.info',
        '/home/liaoth/project/170602_XK/server_result/result/coverage_info/dna_repair/XK-8T-2_cov_extracted.info'
    ]:
        #format_cov_depth_heatmap(path)
        name_id = re.findall('(XK-[28][WT]-?2?)',path)[0]
        avg_depths = normalize_df[normalize_df.loc[:,'Sample ID']==name_id].loc[:,'avg_depth'].values[0]
        hm_data = format_cov_depth_heatmap(path,sum_sig_genes,avg_depths)

        with open(path.replace('_cov_extracted.info','.hm_data'),'w') as f1:
            hm_data.to_csv(f1)

        #plt.figure(figsize=[150,150],dpi=300)
        plt.subplot(1, 6, count)
        if count in [0,2,4]:
            ax = sns.heatmap(hm_data.fillna(0),square=True,xticklabels=5,vmin=0,vmax=800,yticklabels=False)
        else:
            ax = sns.heatmap(hm_data.fillna(0),square=True,xticklabels=5,vmin=0,vmax=800)
        plt.yticks(rotation=0)
        plt.show()
        title_name = ['XK-2W','XK-2T','XK-2T-2','XK-8W','XK-8T','XK-8T-2'][count-1]
        plt.xlabel('%s coverage infomation' % title_name)
        count +=1


    ####using plotly to formatting the data to draw the pictures
    import plotly
    from plotly.graph_objs import *

    fig = plotly.tools.make_subplots(1, 2, shared_yaxes=True, subplot_titles=('XK-27T', 'XK-27W'))
    fig.update(dict(layout=dict(height=1500, width=1500, title='significance genes in Lung Cancer Coverage Heatmap',
                                titlefont={'size': 30})))
    for idx, path in enumerate(['/home/liaoth/project/170801_XK/result/heatmap/XK-27T_sorted_cov.info_all_r_genes',
                                '/home/liaoth/project/170801_XK/result/heatmap/XK-27W_sorted_cov.info_all_r_genes'
                                ]):
        hm_data = df.from_csv(path.replace('_sorted_cov.info_all_r_genes', '.hm_data'))
        fig.data += [Heatmap(z=hm_data.values.tolist(), y=hm_data.index, xaxis='x' + str(idx + 1), zmin=0, zmax=5)]
    plotly.offline.plot(fig)
