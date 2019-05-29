"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
import itertools

import pandas as pd


class fileparser():
    def __init__(self, filename):
        self.filename = filename
        self.df = pd.read_csv(filename, index_col=None)

        self.cols = validate_df(self.df)
        self.df.fillna('')
        self.df = self.df.set_index("sample_name")

    def get_attr(self, col):
        if col not in self.cols:
            raise Exception("attr %s not in input df" % col)
        else:
            return self.df[col].to_dict()

    def is_somatic(self):
        if len(self.df["Somatic"].unique()) > 1:
            return True
        else:
            return False

    def germline_pair(self):
        return self.df.to_dict(orient='index')

    def somatic_pair(self):
        pair_dict = {}
        if self.is_somatic():
            for sample_name in set(self.df.index):
                sub_df = self.df.loc[sample_name, :]
                normal = sub_df.loc[sub_df["Somatic"] == "N"]
                tumor = sub_df.loc[sub_df["Somatic"] == "T"]
                product_combinations = list(itertools.product(range(normal.shape[0]),
                                                              range(tumor.shape[0])))
                if len(product_combinations) > 1:
                    for idx, comb in enumerate(product_combinations):
                        key = "%s-%s" % (sample_name, idx + 1)
                        pair_dict[key] = {}
                        pair_dict[key]["Normal"] = normal.iloc[comb[0], :].to_dict()
                        pair_dict[key]["Tumor"] = tumor.iloc[comb[1], :].to_dict()
                else:
                    comb = product_combinations[0]
                    pair_dict[sample_name] = {}
                    pair_dict[sample_name]["Normal"] = normal.iloc[comb[0], :].to_dict()
                    pair_dict[sample_name]["Tumor"] = tumor.iloc[comb[1], :].to_dict()
            return pair_dict
        else:
            raise Exception("No somatic/pair samples detected")


def validate_df(df):
    template_file = os.path.join(os.path.dirname(__file__), "data_input.template")
    columns_values = open(template_file).read().strip('\n').split(',')

    if set(df.columns) != set(columns_values):
        raise Exception("INPUT file has unknown header")

    # if df["sample_name"].duplicated().any():
    #     raise Exception("sample_name has duplicated.")
    return columns_values


# def pfn(filename, wanted):
#
#     filename = os.path.basename(filename)
#     storged_dict = {}
#     if not PROJECT_NAME:
#         if '_' in filename and '-' in filename:
#             storged_dict['project_name'] = filename[:min(filename.index('_'), filename.index('-'))]
#         else:
#             try:
#                 storged_dict['project_name'] = filename[:filename.index('_')]
#             except:
#                 storged_dict['project_name'] = filename[:filename.index('-')]
#     else:
#         storged_dict['project_name'] = PROJECT_NAME
#     try:
#         storged_dict['sample_name'] = re.findall(sample_name_pattern,
#                                                  filename)[0]
#         storged_dict['sample_ID'] = re.findall(sample_ID_pattern,
#                                                filename)[0]
#     except:
#         print('wrong pattern, please check it.\nOri file name: %s which could not find any given pattern' % (filename))
#
#     # get pair ID (which mean PE1 or PE2)
#     try:
#         storged_dict['pair_ID'] = re.findall(pair_ID_pattern, filename)[0]
#     except:
#         print('wrong pattern, please check it.\nOri file name: %s which could not find any given pattern' % (filename))
#
#     if input_pair_tuple:
#         # if manually input,use it directly.
#         # for some you can't strip simple N/T out to form pair result prefix file.
#         for _pair in input_pair_tuple:
#             # order is important, first one is tumor sample, next is normal sample
#             if _pair[0] in filename:
#                 storged_dict['mt2_for'] = TUMOR_SIG
#                 storged_dict['pair_name'] = storged_dict['project_name'] + '_%s' % input_pair_tuple.index(_pair)
#             elif _pair[1] in filename:
#                 storged_dict['mt2_for'] = NORMAL_SIG
#                 storged_dict['pair_name'] = storged_dict['project_name'] + '_%s' % input_pair_tuple.index(_pair)
#             else:
#                 continue
#     else:
#         # else find sign yes/not
#         if (NORMAL_SIG and TUMOR_SIG):
#             # if both exists
#             try:
#                 storged_dict['mt2_for'] = re.findall(mt2_for_pattern, filename)[0]
#                 storged_dict['pair_name'] = re.findall(pair_name_pattern,filename)[0]
#             except:
#                 print('wrong pattern, please check it.\nOri file name: %s which could not find any given pattern' % (filename))
#
#     if wanted == 'all':
#         return storged_dict
#     else:
#         return storged_dict[wanted]

from main import *
