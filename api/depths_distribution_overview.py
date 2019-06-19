#################################################################################
#### For accessment of Somatic samples
####
#### validated at 190619
#################################################################################
import os
from glob import glob

import click
import pandas as pd
import plotly
import plotly.graph_objs as go

@click.command()
@click.option("--normal", "-n", help="The path with wildcard or not to indicated coverage summary info of Normal samples. End with cov_summary.info")
@click.option("--tumor", "-t", help="The path with wildcard or not to indicated coverage summary info of Tumor samples. End with cov_summary.info")
@click.option("--output", "-o", help="The ouput file name of html")
def main(normal, tumor, output):
    if '*' in normal:
        n_files = list(sorted(glob(normal)))
        t_files = list(sorted(glob(tumor)))
    else:
        n_files = [normal]
        t_files = [tumor]
    # formatted
    n_files = [os.path.abspath(_) for _ in n_files]
    t_files = [os.path.abspath(_) for _ in t_files]

    fig = go.Figure()
    for n_file, t_file in zip(n_files, t_files):
        n_id = os.path.basename(n_file).split('.')[0]
        t_id = os.path.basename(t_file).split('.')[0]
        n_df = pd.read_csv(n_file, sep=',', index_col=0)
        t_df = pd.read_csv(t_file, sep=',', index_col=0)
        fig.add_scatter(x=n_df.loc[:, "depths"],
                        y=n_df.loc[:, "number of pos larger than this depth"],
                        mode='markers+lines',
                        name=n_id,
                        marker=dict(symbol='circle'),
                        legendgroup="normal"
                        )
        fig.add_scatter(x=t_df.loc[:, "depths"],
                        y=t_df.loc[:, "number of pos larger than this depth"],
                        mode='markers+lines',
                        name=t_id,
                        marker=dict(symbol='star-square'),
                        legendgroup="tumor"
                        )
    plotly.offline.plot(fig, filename=output, auto_open=False)


if __name__ == '__main__':
    main()
