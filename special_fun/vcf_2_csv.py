import multiprocessing as mp

import pandas as pd
from tqdm import tqdm


def loader(file_stream, header):
    for line in file_stream:
        yield (_process, (line, header))


def _process(line, header):
    line = line.strip('\n')
    list_line = line.split("\t")
    dict_line = dict(zip(header, list_line))

    INFO_list = dict_line.pop("INFO").split(';')
    INFO_list = [_ for _ in INFO_list if '=' in _]
    INFO_dict = dict([("INFO_" + _.split('=')[0], _.split('=')[1])
                      for _ in INFO_list])
    dict_line.update(INFO_dict)
    samples = header[header.index("FORMAT") + 1:]

    for sample in samples:
        FORMAT_header = ['_'.join([sample, _])
                         for _ in dict_line["FORMAT"].split(":")]
        sample_list = dict_line.pop(sample).split(':')
        sample_dict = dict(zip(FORMAT_header,
                               sample_list))
        dict_line.update(sample_dict)
    dict_line.pop("FORMAT")
    dict_line = {0: dict_line}
    df = pd.DataFrame.from_dict(dict_line, orient='index')
    # df_collect.append(df)
    return df


def exec_func(args):
    func, param = args
    return func(*param)


def vcf2csv(input_vcf, output_tsv, thread=5):
    """

    :param input_vcf:
    :param output_csv:
    :return:
    """
    df_collect = []
    header = []

    if thread == 0 or thread == -1:
        thread = mp.cpu_count()

    with open(input_vcf, 'r') as f1:
        for line in tqdm(f1):
            if line.startswith('##'):
                continue
            elif line.startswith("#"):
                header = line.strip('#').strip('\n').split('\t')
                continue
            break
        df_sorted_cols = header[:header.index("INFO")]
        params = loader(f1, header)
        with mp.Pool(processes=thread) as tp:
            map_r = tp.imap(exec_func, params)
            for df in tqdm(map_r):
                df_collect.append(df)

    final_df = pd.concat(df_collect, axis=0)
    df_sorted_cols += list(sorted(list(final_df.columns.difference(set(df_sorted_cols)))))
    final_df = final_df.reindex(columns=df_sorted_cols)
    final_df.to_csv(output_tsv, index=0, sep='\t')
