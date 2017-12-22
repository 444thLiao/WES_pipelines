import seaborn as sns
from pandas import DataFrame as df
import pandas


query_df = pandas.read_excel('/home/liaoth/project/170602_XK/server_result/result/archived_v2/mRNA/XK_all_rpkm_count_flux_capacitor.xlsx')
all_coding_genes = list(query_df[query_df.gene_type == 'protein_coding'].gene_name)


a= 'XK-8'
mrna_file1 = '/home/liaoth/project/170602_XK/server_result/result/archived_v2/mRNA/degust_%s_two.tsv' % a
mrna_file1 = mrna_file1.replace('XK-8','XK_8')
info_df= df.from_csv(mrna_file1,sep='\t',index_col=False)
new_df = info_df[info_df.gene_name.isin(gene_151)].iloc[:,[0,5,7]]
new_df.index=new_df.gene_name
new_df = new_df.iloc[:,[1,2]]
bucket = []
for gene in list(new_df.index):
    temp_df = df(columns=['genes','sample','values'])
    temp_df.loc[0,:] = [gene,a,new_df.loc[gene,'%sa' % a]]
    temp_df.loc[1, :] = [gene, a+'-2', new_df.loc[gene, '%s-2a' % a]]
    bucket.append(temp_df)
draw_data = pandas.concat(bucket)
sns.pointplot(x='sample',hue='genes',y='values',data=draw_data)
sns.boxplot(x='sample',y='values',data=draw_data)



mrna_file1 = '/home/liaoth/project/170602_XK/server_result/result/archived_v2/mRNA/degust_%s_two.tsv' % a
mrna_file1 = mrna_file1.replace('XK-8','XK_8')
info_df= df.from_csv(mrna_file1,sep='\t',index_col=False)
biggest_10_genes = list(info_df[info_df.gene_name.isin(all_coding_genes)].sort_values('%s-2' % a).gene_name)[-10:]
smallest_10_genes = list(info_df[info_df.gene_name.isin(all_coding_genes)].sort_values('%s-2' % a).gene_name)[:10]

sum_genes = biggest_10_genes+smallest_10_genes