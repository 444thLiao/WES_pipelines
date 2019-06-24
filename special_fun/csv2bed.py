import pandas as pd
import sys
def csv2bed(csv_path, output):
    cache = pd.read_csv(csv_path, index_col=None)
    cache.Start = cache.Start - 1

    sorter = ['chrM',
              'chr1',
              'chr2',
              'chr3',
              'chr4',
              'chr5',
              'chr6',
              'chr7',
              'chr8',
              'chr9',
              'chr10',
              'chr11',
              'chr12',
              'chr13',
              'chr14',
              'chr15',
              'chr16',
              'chr17',
              'chr18',
              'chr19',
              'chr20',
              'chr21',
              'chr22',
              'chrX',
              'chrY',
              'chr1_gl000191_random',
              'chr1_gl000192_random',
              'chr4_ctg9_hap1',
              'chr4_gl000193_random',
              'chr4_gl000194_random',
              'chr6_apd_hap1',
              'chr6_cox_hap2',
              'chr6_dbb_hap3',
              'chr6_mann_hap4',
              'chr6_mcf_hap5',
              'chr6_qbl_hap6',
              'chr6_ssto_hap7',
              'chr7_gl000195_random',
              'chr8_gl000196_random',
              'chr8_gl000197_random',
              'chr9_gl000198_random',
              'chr9_gl000199_random',
              'chr9_gl000200_random',
              'chr9_gl000201_random',
              'chr11_gl000202_random',
              'chr17_ctg5_hap1',
              'chr17_gl000203_random',
              'chr17_gl000204_random',
              'chr17_gl000205_random',
              'chr17_gl000206_random',
              'chr18_gl000207_random',
              'chr19_gl000208_random',
              'chr19_gl000209_random',
              'chr21_gl000210_random',
              'chrUn_gl000211',
              'chrUn_gl000212',
              'chrUn_gl000213',
              'chrUn_gl000214',
              'chrUn_gl000215',
              'chrUn_gl000216',
              'chrUn_gl000217',
              'chrUn_gl000218',
              'chrUn_gl000219',
              'chrUn_gl000220',
              'chrUn_gl000221',
              'chrUn_gl000222',
              'chrUn_gl000223',
              'chrUn_gl000224',
              'chrUn_gl000225',
              'chrUn_gl000226',
              'chrUn_gl000227',
              'chrUn_gl000228',
              'chrUn_gl000229',
              'chrUn_gl000230',
              'chrUn_gl000231',
              'chrUn_gl000232',
              'chrUn_gl000233',
              'chrUn_gl000234',
              'chrUn_gl000235',
              'chrUn_gl000236',
              'chrUn_gl000237',
              'chrUn_gl000238',
              'chrUn_gl000239',
              'chrUn_gl000240',
              'chrUn_gl000241',
              'chrUn_gl000242',
              'chrUn_gl000243',
              'chrUn_gl000244',
              'chrUn_gl000245',
              'chrUn_gl000246',
              'chrUn_gl000247',
              'chrUn_gl000248',
              'chrUn_gl000249']
    cache.Chr = cache.Chr.astype('category')
    cache.Chr.cat.set_categories(sorter, inplace=True)

    cache.sort_values(['Chr', 'Start', 'End'], inplace=True)

    with open(output, 'w') as f1:
        cache.iloc[:, :3].to_csv(f1, index=False, header=None, sep='\t')
    return output

if __name__ == '__main__':
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    csv2bed(input_path, output_path)
