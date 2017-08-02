from pandas import DataFrame as df

def join_2_state_df(left_csv,right_csv):
    left_one = df.from_csv(left_csv,index_col=False)
    right_one = df.from_csv(right_csv,index_col=False)
    left_one.index = [';'.join([str(_i) for _i in list(left_one.iloc[_idx, :5])]) for _idx in list(left_one.index)]
    right_one.index = [';'.join([str(_i) for _i in list(right_one.iloc[_idx, :5])]) for _idx in list(right_one.index)]
    venn2_unweighted([set(left_one.index), set(right_one.index)], set_labels=['XK-2', 'XK-2-2'])