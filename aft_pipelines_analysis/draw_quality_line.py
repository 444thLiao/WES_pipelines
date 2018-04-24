import plotly
from plotly.graph_objs import *
import glob,tqdm,os
import pandas as pd
import sys

def draw_coverage_depths(dir_info,NORMAL_SIG,TUMOR_SIG):
    draw_data = []
    layout = dict(yaxis=dict(exponentformat='e'))
    all_info = list(glob.glob(dir_info+'/*cov_summary.info'))
    for data_path in all_info:
        data_df = pd.read_csv(data_path,sep='\t',index_col=0)
        samples_names = os.path.basename(data_path).split('cov_summary')[0].strip('_')
        if NORMAL_SIG in samples_names:
            draw_data.append(
                Scatter(x=data_df.loc[:, 'depths'], y=data_df.loc[:, 'coverage'], mode='markers+lines',
                        name=samples_names,
                        marker=dict(symbol='circle'),
                        legendgroup=samples_names.replace(NORMAL_SIG, '').replace(TUMOR_SIG, '')))
        elif TUMOR_SIG in samples_names:
            draw_data.append(
                Scatter(x=data_df.loc[:,'depths'],y=data_df.loc[:,'coverage'],mode='markers+lines',
                        name=samples_names,
                        marker=dict(symbol='star-square'),
                        legendgroup=samples_names.replace(NORMAL_SIG,'').replace('%s-2' % TUMOR_SIG, '').replace('T','')))
    plotly.offline.plot(dict(data=draw_data, layout=layout), filename=dir_info + '/depth_dis_lines.html')


if __name__ == '__main__':
    if len(sys.argv) >=2:
        project_name = sys.argv[-1]
    else:
        project_name = '180321'
    input_path = '/home/liaoth/data2/project/XK_WES/%s/output/info_summary/' % project_name
    draw_data = []
    layout=dict(yaxis=dict(exponentformat='e'))
    for data in tqdm.tqdm(glob.glob(input_path+'/*cov_summary.info')):
        data_df = pd.read_csv(data,sep='\t',index_col=0)
        samples_names = os.path.basename(data).split('cov_summary')[0].strip('_')
        if NORMAL_SIG in samples_names:
            draw_data.append(
                Scatter(x=data_df.loc[:, 'depths'], y=data_df.loc[:, 'coverage'], mode='markers+lines',
                        name=samples_names,
                        marker=dict(symbol='circle'),
                        legendgroup=samples_names.replace('W', '').replace('T', '')))
        elif TUMOR_SIG in samples_names:
            draw_data.append(
                Scatter(x=data_df.loc[:,'depths'],y=data_df.loc[:,'coverage'],mode='markers+lines',
                        name=samples_names,
                        marker=dict(symbol='star-square'),
                        legendgroup=samples_names.replace('W','').replace('T-2', '').replace('T','')))

    plotly.offline.plot(dict(data=draw_data,layout=layout),filename=input_path+'/depth_dis_lines.html')