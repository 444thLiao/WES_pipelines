
from pandas import DataFrame as df
from collections import Counter
import tqdm
def merge_overlap_tuple(alist):
    """
    input a bit list with lots of tuple inside it.
    tuple have pos.like  (413, 414),
 (415, 416),
 (418, 419),
 (418, 419, 420),
 (419, 420, 421),
 (420, 421, 422),
 (421, 422),
 (423, 424)
    must be in order.
    we could merge  (418, 419),(418, 419, 420),(419, 420, 421),(420, 421, 422),(421, 422) into one (418,419,420,421,422)
    :param alist:
    :return:
    """
    new_bucket = []
    past_t = alist[0]
    for current_t in alist[1:]:
        if set(current_t).intersection(set(past_t)):
            new_bucket = new_bucket[:-1]
            new_bucket.append(tuple(sorted(set(list(current_t)+list(past_t)))))
        else:
            new_bucket.append(current_t)
        past_t = current_t
    return new_bucket


def detect_possible(csv_path,length=150/2):
    bucket = []
    cache_df = df.from_csv(csv_path,index_col=False)
    cache_df.sort_values(['Chr','Start'],inplace=True)
    cache_df.index = range(len(cache_df))
    chrs = sorted(list(set(cache_df.Chr)))
    for _chr in chrs:
        sub_cache = cache_df[cache_df.Chr == _chr]
        for _idx in range(len(sub_cache)):
            query_pos = int(sub_cache.iloc[_idx,1])

            right_one = sub_cache[(sub_cache.Start > query_pos-length) & (sub_cache.Start < query_pos+length)]
            if len(right_one) > 1 and tuple(right_one.index) not in bucket:
                bucket.append(tuple(sorted(right_one.index)))
    bucket = merge_overlap_tuple(bucket)
    #longer and merge neighbour snv.
    return bucket

def parse_startend(list_of_pos):
    starts = []
    ends = []
    for i in list_of_pos:
        starts.append(i[1])
        ends.append(i[2])
    return (min(starts+ends),max(starts+ends))

import pysam
from collections import defaultdict
def fetch_possible_reads(bam,list_of_pos):
    """

    :param bam:
    :param list_of_pos: a list with tuple. [('chr',start,end,ref,alt),()]
    :return:
    """
    bam_file = pysam.AlignmentFile(bam)
    start,end = parse_startend(list_of_pos)
    reads = bam_file.fetch(list_of_pos[0][0],int(start),int(end))
    reads = [_i for _i in reads]
    snv2reads = defaultdict(list)
    for snv in list_of_pos:
        if snv[3].upper() in ['A', 'C', 'T', 'G'] and snv[4].upper() in ['A', 'C', 'T', 'G']:
            for read in reads:
                if read.is_unmapped:
                    continue
                if snv[1] in range(read.aend-read.alen,read.aend):
                    try:
                        _idx = read.get_reference_positions().index(snv[1])-1
                        #if this pos is a deletion, it wiil disapper from this function result.
                    except:
                        continue
                    fetch_alt = read.query_alignment_sequence[_idx]
                    fetch_ref = read.get_reference_sequence()[_idx]
                    if fetch_ref.upper() == snv[3] and fetch_alt.upper() == snv[4]:
                        snv2reads[snv].append(read)

        else:
            continue
    return snv2reads

def cal_cross_both_pos(snv1,snv2,reads):
    bucket = []
    for read in reads:
        if snv1[1] in range(read.aend - read.alen, read.aend) and snv2[1] in range(read.aend - read.alen, read.aend):
            bucket.append(read)
    return len(bucket)

def parse_possible_reads(snv2reads):
    pair_compare = {}
    for c_k,c_v in snv2reads.items():
        for _k,_v in snv2reads.items():
            if _k != c_k and tuple(sorted([c_k,_k])) not in pair_compare.keys():
                count = len(set(c_v).intersection(set(_v)))
                pair_compare[tuple(sorted([c_k,_k]))] = {'same reads':count,
                                                         'each reads':{c_k:len(c_v),_k:len(_v)},
                                                         'cross both pos':{c_k:cal_cross_both_pos(c_k,_k,c_v),
                                                                           _k:cal_cross_both_pos(c_k,_k,_v)}}
    return pair_compare

def construct_haplotype(matrix_of_haplotype):
    """
    It is one-one neightbour compare, so we will lost some far-distance haplotype.
    :param matrix_of_haplotype:
    :return:
    """
    bucket = []
    for _k in matrix_of_haplotype.keys():
        bucket += list(_k)
    bucket = sorted(list(set(bucket)))
    ### already in ordered list.
    haplotype = []
    for _idx,snv in enumerate(bucket):
        if _idx != len(bucket)-1:
            current_relation = matrix_of_haplotype[tuple(sorted([snv,bucket[_idx+1]]))]
            if current_relation['cross both pos'].values()[0] == 0 or current_relation['cross both pos'].values()[1] == 0:
                continue
            if current_relation['same reads']/float(current_relation['cross both pos'].values()[0]) > 0.8 or \
                current_relation['same reads'] / float(current_relation['cross both pos'].values()[1]) > 0.8:
                try:
                    if haplotype[-1] != snv:
                        haplotype += ['|',snv,bucket[_idx+1]]
                    else:
                        haplotype += [bucket[_idx+1]]
                except:
                    haplotype += [snv, bucket[_idx + 1]]
            else:
                try:
                    if haplotype[-1] != snv:
                        haplotype += ['|',snv,'|',bucket[_idx+1]]
                    else:
                        haplotype += ['|',bucket[_idx+1]]
                except:
                    haplotype += [snv, '|', bucket[_idx + 1]]
    return haplotype

def construct_haplotype_v2(matrix_of_haplotype):
    bucket = []
    for _k in matrix_of_haplotype.keys():
        bucket += list(_k)
    bucket = sorted(list(set(bucket)))
    ### already in ordered list.
    pre_haplotype = []
    haplotype = []
    for _idx,snv in enumerate(bucket[:-1]):
        #in case take last one.
        _cache = [snv]
        for new_snv in bucket[_idx+1:]:
            current_relation = matrix_of_haplotype[tuple(sorted([snv,new_snv]))]
            if current_relation['cross both pos'].values()[0] == 0 or current_relation['cross both pos'].values()[1] == 0:
                continue
            if current_relation['same reads']/float(current_relation['cross both pos'].values()[0]) > 0.8 or \
                current_relation['same reads'] / float(current_relation['cross both pos'].values()[1]) > 0.8:
                _cache += [new_snv]
        pre_haplotype.append(_cache)
    try:
        pre_haplotype.append([bucket[-1]])
    except:
        pass

    while True:
        if len(pre_haplotype) >= 2:
            first = pre_haplotype[0]
            for each in pre_haplotype[1:]:
                if set(first).intersection(set(each)):
                    pre_haplotype[0] = list(set(first + each))
                    pre_haplotype.remove(each)
                    break
                if each == pre_haplotype[-1]:
                    pre_haplotype = pre_haplotype[1:]+[pre_haplotype[0]]

        count_for = []
        for i in pre_haplotype:
            count_for+=i
        count_for = Counter(count_for)
        if sum(count_for.values()) == len(bucket):
            break
    for _i in pre_haplotype:
        haplotype+=_i
        haplotype.append('|')
    haplotype = haplotype[:-1]
    return haplotype




def find_haplotype(csv,bam):
    _list = detect_possible(csv)
    _df = df.from_csv(csv,index_col=False)
    _df.sort_values(['Chr','Start'],inplace=True)
    _df.index = range(len(_df))
    bucket = []
    for each in tqdm.tqdm(_list):
        _cache = [tuple(_df.iloc[_idx,:5])for _idx in each]
        snv2reads = fetch_possible_reads(bam,_cache)
        matrix_of_haplotype = parse_possible_reads(snv2reads)
        haplotype_summary = construct_haplotype_v2(matrix_of_haplotype)

        bucket.append(haplotype_summary)
        bucket = [_i for _i in bucket if _i]
    return bucket
#
# fetch_possible_reads('/home/liaoth/project_formal/170602_XK/output/XK_result/XK-8T-2/XK-8T-2.recal_reads.bam',[('chr1', 16386447, 16386447, 'G', 'C'),('chr1', 16386495, 16386495, 'C', 'T')])


def haplotype_report(haplotype_bucket,output_file):
    report_str = 'Chr\tStart\tEnd\tRef\tAlt\tGroup_ID\tHaplotype_ID\n'
    count_group = 0

    for SNVs in tqdm.tqdm(haplotype_bucket):
        count_group +=1
        haplotype_id = 0
        for SNV in SNVs:
            if SNV != '|':
                report_str+='\t'.join([str(_i) for _i in SNV])
                report_str+='\t'+ str(count_group)+'\t'+str(haplotype_id)+'\n'
            else:
                haplotype_id+=1
    with open(output_file,'w') as f1:
        f1.write(report_str)


if __name__ == '__main__':
    import sys
    csv_file = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    haplotype_report(find_haplotype(csv_file, bam_file), output_file)
