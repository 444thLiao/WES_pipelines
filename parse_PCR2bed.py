from pandas import DataFrame as df


"""
Function for convert a blastn result(e.g PCR gene blastn result) to a formatted bed file.

:type Formatted-for function.
"""

def parse_PCR2BED(in_path,out_path,total_gene):

    process_one = df.from_csv(in_path, sep='\t', header=None, index_col=False)
    process_one.columns = blastn_header

    process_one = process_one.sort_values('pident', ascending=False)
    process_one.index = range(len(process_one))

    total_gene = int(total_gene)
    if len(process_one[process_one.pident == 100]) == total_gene:
        result = process_one.iloc[:total_gene, [0, 1, 3, 8, 9]]
        result_txt = bed_file_fmt
        for nrow in range(len(result)):
            result_txt += result.iloc[nrow, 1] + '\t'
            if result.iloc[nrow, 3] < result.iloc[nrow, 4]:
                strand = '+'
            else:
                strand = '-'
            start = str(min(result.iloc[nrow, 3], result.iloc[nrow, 4]))
            end = str(max(result.iloc[nrow, 3], result.iloc[nrow, 4]))
            result_txt += '%s\t%s\t%s\t%s\n' % (start, end, strand, result.iloc[nrow, 0])

        with open(out_path, 'w') as f1:
            f1.write(result_txt)
    else:
        print 'Conflict need to handle in parse_PCR2BED.'

from main import *