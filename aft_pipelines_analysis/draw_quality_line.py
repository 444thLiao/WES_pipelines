import plotly
from plotly.graph_objs import *
import glob,tqdm,os
import pandas as pd
input_path = '/home/liaoth/project/180104_XK/gpz_server/quality_assessment/'
draw_data = []
layout=dict(yaxis=dict(exponentformat='e'))
for data in tqdm.tqdm(glob.glob(input_path+'/*cov_summary.info')):
    data_df = pd.read_csv(data,sep='\t',index_col=0)
    samples_names = os.path.basename(data).split('cov_summary')[0].strip('_')
    draw_data.append(Scatter(x=data_df.loc[:,'depths'],y=data_df.loc[:,'coverage'],mode='markers+lines',name=samples_names))

plotly.offline.plot(dict(data=draw_data,layout=layout),filename=input_path+'/depth_dis_lines.html')