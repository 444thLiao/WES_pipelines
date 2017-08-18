#from __future__ import absolute_import
import pandas
def extract2dict(Otherinfo):
    _info = Otherinfo.split('\t')
    _n = _info[-2].split(':')
    _v = _info[-1].split(':')
    _query = dict(zip(_n, _v))
    return _query

def count_indel(file_df):
    count = []
    for _index in list(file_df.index):
        if file_df.loc[_index,'Ref'] in ['A','C','T','G'] and file_df.loc[_index,'Alt'] in ['A','C','T','G']:
            count.append('SNP')
        elif file_df.loc[_index,'Ref'] == '-' and len(file_df.loc[_index,'Alt']) >=1:
            count.append('insert')
        elif file_df.loc[_index,'Alt'] == '-' and file_df.loc[_index,'Ref'] >=1:
            count.append('del')
        else:
            print file_df.loc[_index,'Alt'],file_df.loc[_index,'Ref']
    result = {}
    for _i in set(count):
        result[_i] = count.count(_i)
    result['Indel'] = result['del']+result['insert']
    return result


def counter_pass(file_df):
    from collections import Counter
    list_pass = []
    for _index in list(file_df.index):
        list_pass+=file_df.loc[_index,'Otherinfo'].split('\t')[6].split(';')

    return Counter(list_pass)


def construct_pos_list(bed_df,coord = 0):
    '''
    :param bed_df: a pandas.DataFrame as input, it need to have #chrom,start,end three base columns.
    :param cor: coordiante you use, basically bed,bam,bcf is using 0-coordinate.
    :return:
    '''
    from collections import defaultdict
    pos_info_dict = defaultdict(list)

    for _index in list(bed_df.index):
        chr = bed_df.loc[_index, '#chrom']
        start = bed_df.loc[_index, 'start']
        end = bed_df.loc[_index, 'end']
        if coord == 0:
            intervals = range(start+1, end + 1) # turn to 1-corordiante.
        elif coord == 1:
            intervals = range(start,)
        else:
            #print 'cor need to be 0/1, meaning the corordinate you use.'
            raise SyntaxError,'cor need to be 0/1, meaning the corordinate you use.'
        pos_info_dict[chr] += intervals
    return pos_info_dict


from pandas import DataFrame as df
import pysam
from Whole_pipelines.setting import REF_file_path
def cal_fun(bam_path, bed_file):
    '''
    Using a input bed file and a bam file, cal each pos cov_info from the intervals in bed.

    :param bam_path:
    :param bed_file:
    :return: writing the cov_info to file besides the bam file, with 1-coordinate.
    '''
    bed_file = df.from_csv(bed_file, index_col=False, sep='\t')
    try:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
    except:
        print "bamefile need to cal doesn't exist"
        raise IOError

    fastafile = pysam.FastaFile(filename=REF_file_path)
    ####define each file path.
    result = 'Gene\tChr\tPosistion\tReference\tbase\tA\tC\tG\tT\tA_Rate\tC_Rate\tG_Rate\tT_Rate\t1\t2\t3\t4\n'
    #define columns.
    for i in range(len(bed_file)):
        # iterator
        Chr = bed_file.iloc[i, 0]
        start = min(int(bed_file.iloc[i, 2]), int(bed_file.iloc[i, 1]))
        end = max(int(bed_file.iloc[i, 2]), int(bed_file.iloc[i, 1]))
        Gene_name = str(bed_file.iloc[i, 3])
        #fetch basic info.
        coverage_ACGT = bamfile.count_coverage(Chr, start - 1, end, read_callback='nofilter',
                                               quality_threshold=0) # from 1- to 0-
        base_counter = dict(zip(['A', 'C', 'G', 'T'], coverage_ACGT))
        # make it a dict according preset order.
        ref_base_str = fastafile.fetch(Chr, start - 1, end) #need to input 0- start/end.

        if sum(sum(j) for j in coverage_ACGT) != 0:
            # if this intervals doesn't have any reads, we ignore it.
            for base_n in range(start, end + 1):
                n_read = int(sum([k[range(start, end + 1).index(base_n)] for k in coverage_ACGT]))
                # total base num
                if n_read == 0:
                    continue
                    #if position here didn't have any base, then pass this pos.
                result += Gene_name + '\t' + Chr + '\t' + str(base_n) + '\t' + ref_base_str[
                    range(start, end + 1).index(base_n)] + '\t' + '\t'
                for base in ['A', 'C', 'G', 'T']:
                    result += str(base_counter[base][range(start, end + 1).index(base_n)]) + '\t'
                    #write the A/T/G/C num.

                for base in ['A', 'C', 'G', 'T']:
                    result += str(
                            round(float(base_counter[base][range(start, end + 1).index(base_n)]) / n_read, 4)) + '\t'
                    # wirte the rate of A/T/C/G RATE
                little_rank_list = [(each_one[1][range(start, end + 1).index(base_n)], each_one[0]) for each_one in
                                  base_counter.items()]
                #construct a list to sort A/T/C/G each base num.
                result += '\t'.join([ranked_base[1] for ranked_base in sorted(little_rank_list, reverse=True)]) + '\n'
                #make it columns and tag a '\n' sign.
            if result[-1:] != '\n':
                result += '\n'
            #in case the last base didn't have any base, it will continue, so need to check the last sign.
        else:
            pass
    with open(bam_path.partition('.')[0]+'_cov.info', 'w') as f1:
        f1.write(result)
    print 'Cal cov info complete'


def pfn(filename, wanted):
    import setting
    import re
    '''
    Basic function for parse input sample name to formatted format in order to import into pipelines.
    :type Parse_Function

    :param filename:
    :param wanted:
    :return:
    '''
    if '/' in filename:
        filename = filename.rpartition('/')[2]
    storged_dict = {}
    if '_' in filename and '-' in filename:
        storged_dict['project_name'] = filename[:min(filename.index('_'), filename.index('-'))]
    else:
        try:
            storged_dict['project_name'] = filename[:filename.index('_')]
        except:
            storged_dict['project_name'] = filename[:filename.index('-')]

    #storged_dict['sample_name'] = filename.partition('_S')[0]
    storged_dict['sample_name'] = re.findall(setting.sample_name_pattern,filename)[0]

    new_filename = filename[filename.index('_') + 1:]
    storged_dict['sample_ID'] = re.findall(setting.sample_ID_pattern, new_filename)[0]
    if 'R' in new_filename:
        try:
            storged_dict['pair_ID'] = re.findall(setting.pair_ID_pattern, new_filename)[0]
        except:
            print new_filename
    if setting.NORMAL_SIG in filename or setting.TUMOR_SIG in filename:
        storged_dict['mt2_for'] = re.findall(setting.mt2_for_pattern, filename)[0]
        storged_dict['pair_name'] = re.findall(setting.pair_name_pattern,filename)[0]

    if wanted == 'all':
        return storged_dict
    else:
        return storged_dict[wanted]

def csv_compare(csv1,csv2,compare_type = 'diff'):
    info_df1 = df.from_csv(csv1,sep=',',index_col=False)
    info_df2 = df.from_csv(csv2,sep=',',index_col=False)
    info_df1_index = [';'.join([str(_i) for _i in list(info_df1.iloc[_idx, :5])]) for _idx in list(info_df1.index)]
    info_df2_index = [';'.join([str(_i) for _i in list(info_df2.iloc[_idx, :5])]) for _idx in list(info_df2.index)]

    info_df1.index = info_df1_index
    info_df2.index = info_df2_index

    if compare_type=='diff':
        return {'csv1':info_df1.ix[set(info_df1_index).difference(set(info_df2_index))],
                'csv2':info_df2.ix[set(info_df2_index).difference(set(info_df1_index))]}
    elif compare_type=='add':
        new_df = pandas.concat([info_df1, info_df2])
        new_df.index= range(len(new_df))
        return new_df
    #TODO

def typeing_with_gene(extracted_info_df,query_df,query_fields='Gene Name',fetch_field=['Mechanism']):
    gene_list = extracted_info_df.loc[:,'Gene.refGene']
    index_list = list(extracted_info_df.index)
    gene_list = [_i.split(';') for _i in gene_list]

    for idx,genes in zip(index_list,gene_list):
        extracted_mechanism = query_df[query_df.loc[:,query_fields].isin(genes)].loc[:,fetch_field].values.tolist()
        extracted_mechanism = list(set([_m[0] for _m in extracted_mechanism]))
        if 'mechanism' not in extracted_info_df.columns:
            extracted_info_df.loc[:,'mechanism'] = ''

        extracted_info_df.loc[idx,'mechanism'] = ';'.join(extracted_mechanism)
    return extracted_info_df