import plotly
from plotly.graph_objs import *
import glob,tqdm,os
import pandas as pd
import sys

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
    if 'W' in samples_names:
        draw_data.append(
            Scatter(x=data_df.loc[:, 'depths'], y=data_df.loc[:, 'coverage'], mode='markers+lines',
                    name=samples_names,
                    marker=dict(symbol='circle'),
                    legendgroup=samples_names.replace('W', '').replace('T', '')))
    elif 'T' in samples_names:
        draw_data.append(
            Scatter(x=data_df.loc[:,'depths'],y=data_df.loc[:,'coverage'],mode='markers+lines',
                    name=samples_names,
                    marker=dict(symbol='star-square'),
                    legendgroup=samples_names.replace('W','').replace('T-2', '').replace('T','')))

plotly.offline.plot(dict(data=draw_data,layout=layout),filename=input_path+'/depth_dis_lines.html')