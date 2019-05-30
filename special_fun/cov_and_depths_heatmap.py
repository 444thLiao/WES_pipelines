import pandas as pd
from pandas import DataFrame as df
import glob, os, tqdm, re
import plotly
from plotly.graph_objs import *


import pandas

normalize_df = pandas.read_excel('/home/liaoth/project/180104_XK/180104_XK总体评估结果.xlsx')
avg_depth = {s_n: depth for s_n, depth in normalize_df.loc[:, ['Sample ID', 'avg_depth']].values}
print(avg_depth)
avg_depth['25W'] = 178.105

# count = 1
if __name__ == '__main__':
    names = {}
    for path in tqdm.tqdm(glob.glob("/home/liaoth/project/180104_XK/gpz_server/heatmap/*.heatmap")):
        name_id = re.findall('(XK-[0-9]{2}[A-Z])', path)[0]
        heatmap_data = pd.read_csv(path, index_col=0)
        names[name_id] = heatmap_data

    fig = plotly.tools.make_subplots(1, len(names), shared_yaxes=True, subplot_titles=tuple(sorted(names)))
    ####using plotly to formatting the data to draw the pictures


    fig.update(dict(layout=dict(height=1500, width=1500, title='significance genes in Lung Cancer Coverage Heatmap',
                                titlefont={'size': 30})))
    for idx, sample_name in enumerate(sorted(names.keys())):
        hm_data = names[sample_name]
        fig.data += [Heatmap(z=hm_data.values / avg_depth[sample_name.replace('XK-', '')], y=hm_data.index,
                             xaxis='x' + str(idx + 1), zmin=0, zmax=5)]
    plotly.offline.plot(fig, filename='/home/liaoth/project/180104_XK/gpz_server/heatmap/summary_heatmap.html')
