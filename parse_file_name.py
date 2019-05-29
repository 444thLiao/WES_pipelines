"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
import itertools

import pandas as pd
import os

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
        t_dict = self.df.to_dict(orient='index')
        for k,d in t_dict.items():
            d["SampleID"] = k
        return t_dict

    def somatic_pair(self):
        """

        :return: nid + tid: {Normal:{info}},{Tumor:{info}},{nid + tid}
        """
        pair_dict = {}
        if self.is_somatic():
            for source_name in set(self.df.source_name):
                # todo : use group by
                sub_df = self.df.loc[self.df.source_name==source_name, :]
                if len(sub_df.shape) == 1:
                    # only one source_name occur
                    # it may not a tumor-normal pair experiment, pass it
                    continue
                normal = sub_df.loc[sub_df["Somatic"] == "N"]
                tumor = sub_df.loc[sub_df["Somatic"] == "T"]
                product_combinations = list(itertools.product(range(normal.shape[0]),
                                                              range(tumor.shape[0])))
                # it may have multiple pair of samples. e.g. 2 Normal samples vs 1 Tumor samples
                for comb in product_combinations:
                    nid = normal.index[comb[0]]
                    tid = tumor.index[comb[1]]
                    key = "%s + %s" % (nid,
                                       tid)
                    pair_dict[key] = {}
                    pair_dict[key]["Normal"] = normal.iloc[comb[0], :].to_dict()
                    pair_dict[key]["Normal"]["SampleID"] = nid
                    pair_dict[key]["Tumor"] = tumor.iloc[comb[1], :].to_dict()
                    pair_dict[key]["Normal"]["SampleID"] = tid

            return pair_dict
        else:
            raise Exception("No somatic/pair samples detected")


def validate_df(df):
    template_file = os.path.join(os.path.dirname(__file__),
                                 "data_input.template")
    columns_values = open(template_file).read().strip('\n').split(',')

    if set(df.columns) != set(columns_values):
        raise Exception("INPUT file has unknown header")

    if df["sample_name"].duplicated().any():
        raise Exception("sample_name has duplicated.")
    return columns_values
