import re,glob,pandas
#from pandas import DataFrame as df
#import seaborn as sns
import time
def cov_depth(cov_info):
    """
    cal the info from cal_Cov.py and produce a plot
    :param cov_info:
    :return:
    """
    raw_info = pandas.read_csv(cov_info,sep='\t',index_col=False,engine='python')
    raw_info.base = list(raw_info.loc[:, ['A', 'T', 'C', 'G']].sum(1))

    d1 = len(raw_info[raw_info.base > 1])
    d5 = len(raw_info[raw_info.base > 5])
    d10 = len(raw_info[raw_info.base > 10])
    d20 = len(raw_info[raw_info.base > 20])
    result = sum(list(raw_info.base))
    return result,[d1,d5,d10,d20]



df_bucket = []
for path in glob.glob('/home/liaoth/project_formal/170123_XK/output/XK_result/XK-*/XK-*_sorted_sorted_cov.info'):
    print 'processing......',path
    t1=time.time()
    total_bp_in_target,depth_per = cov_depth(path)
    print total_bp_in_target,depth_per
    print 'down...one.....used',time.time()-t1

# a = pandas.concat(df_bucket)
# with open('/home/liaoth/project_formal/170123_XK/cov_for_plot','w') as f1:
#     a.to_csv(f1)

#
# ax = sns.stripplot(x='dep',y='cov',hue='sample',data=a)
#
# depth_name2 = []
# for _idx,_v in enumerate(depths_name):
#     if int(_v.replace('X','')) % 50 ==0:
#         depth_name2.append(_v)
#     else:
#         depth_name2.append('')
# ax.set_xticklabels(depth_name2)
# ax.figure.savefig('/home/liaoth/Desktop/XK_dep_cov_rela.png')