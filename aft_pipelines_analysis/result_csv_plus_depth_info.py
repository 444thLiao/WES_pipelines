import pandas
from pandas import DataFrame as df


"""
Function to add extract column from depth_info file to csv_result file.

Type: Formatted-for function.
"""
def plus_depth_info(depth_info, csv_result):
    csv_f = df.from_csv(csv_result, index_col=False)
    depth_f = df.from_csv(depth_info, index_col=False, sep='\t')

    extract_part_list = []
    otherinfo_df = df(columns=['GT', 'AD', 'AF', 'Filter'])
    for i in range(len(csv_f)):
        chr_query = csv_f.iloc[i, 0]
        Start_query = int(csv_f.iloc[i, 1])
        End_query = int(csv_f.iloc[i, 2])
        cache1 = depth_f[depth_f.loc[:, 'Chr'] == chr_query]

        otherinfo = csv_f.loc[:, 'Otherinfo'][i]

        alt_query = csv_f.iloc[i, 4]
        ref_query = csv_f.iloc[i, 3]
        if alt_query == '-' or ref_query == '-':
            continue

        otherinfo_s = otherinfo.split('\t')

        Filter = otherinfo_s[6]
        GENOTYPE_DICT = dict(zip(otherinfo_s[8].split(':'), otherinfo_s[9].split(':')))
        GT = GENOTYPE_DICT['GT']
        AD = GENOTYPE_DICT['AD']
        try:
            AF = float(GENOTYPE_DICT['AF'])
        except:
            AF = None

        otherinfo_df.loc[i] = [GT, AD, AF, Filter]

        if End_query == Start_query:
            cache2 = cache1[cache1.loc[:, 'Posistion'] == Start_query]
            extract_part = cache2.loc[:, columns_need_to_extract]
            extract_part.index = [i]
            extract_part_list.append(extract_part)
        else:
            pass
            #indel doesn't cal any base info.

    merge_df = pandas.concat(extract_part_list)
    finish_csv1 = csv_f.join(merge_df)
    finish_csv2 = finish_csv1.join(otherinfo_df)

    finish_csv2 = finish_csv2.loc[:, column_retained]
    #reformatted.
    return finish_csv2

from main import *