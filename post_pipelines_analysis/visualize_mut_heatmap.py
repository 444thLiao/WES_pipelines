import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform


# get data
proteins_df = pd.read_csv('/home/liaoth/project/brca_171208/output/mut_counts.csv',index_col=0)
# data_array = data.values
# labels = proteins_df.index.values

figure = ff.create_dendrogram(proteins_df.T.values, orientation='bottom')
for i in range(len(figure['data'])):
    figure['data'][i]['yaxis'] = 'y2'

# Create Side Dendrogram
dendro_side = ff.create_dendrogram(proteins_df.values, orientation='right')
for i in range(len(dendro_side['data'])):
    dendro_side['data'][i]['xaxis'] = 'x2'

# Add Side Dendrogram Data to Figure
figure['data'].extend(dendro_side['data'])

dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
x_order = list(map(int, dendro_leaves))
y_order = list(map(int, figure['layout']['xaxis']['ticktext']))
figure['layout']['xaxis']['ticktext'] = proteins_df.columns.values[y_order]
figure['layout']['yaxis']['ticktext'] = proteins_df.index.values[x_order]
figure['layout']['yaxis']['tickvals'] = dendro_side['layout']['yaxis']['tickvals']
heatmap = go.Heatmap(y=proteins_df.iloc[x_order, y_order].index.values, z=proteins_df.iloc[x_order, y_order].values,
                     colorscale=[[0, '#FFFFFF'], [0.5, '#FFFFFF'],
                                 [0.5, '#bb2a34'], [1, '#bb2a34']],
                     )

heatmap['x'] = figure['layout']['xaxis']['tickvals']
heatmap['y'] = dendro_side['layout']['yaxis']['tickvals']

# Add Heatmap Data to Figure
figure['data'].append(heatmap)

figure['layout'].update({'width':2500, 'height':1500,
                         'showlegend':False, 'hovermode': 'closest',
                         'margin':{'b':120}
                         })
# Edit xaxis
figure['layout']['xaxis'].update({'domain': [.15, 1],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'ticks':""})
# Edit xaxis2
figure['layout'].update({'xaxis2': {'domain': [0, .10],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})

# Edit yaxis
figure['layout']['yaxis'].update({'domain': [0, .85],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
'showticklabels': False,
                                  'ticks': ""})
# Edit yaxis2
figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})
plotly.offline.plot(figure)