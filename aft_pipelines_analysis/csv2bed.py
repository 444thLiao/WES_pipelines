from pandas import DataFrame as df
def csv2bed(csv_path,output):
    cache = df.from_csv(csv_path,index_col=False)
    cache.Start = cache.Start-1
    cache.sort_values(['Chr','Start','End'],inplace=True)
    with open(output,'w') as f1:
        cache.iloc[:,:3].to_csv(f1,index=False,header=None,sep='\t')