"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
import itertools
import os

import pandas as pd

import setting as config


class fileparser():
    def __init__(self, filename):
        filename = os.path.abspath(filename)

        self.df = pd.read_csv(filename, sep='\t', index_col=None, dtype=str)

        self.cols, self.df = validate_df(self.df, filename)
        self.df = self.df.set_index("sample ID")

    def get_attr(self, col):
        if col == self.df.index.name:
            return list(self.df.index)
        if col not in self.cols:
            raise Exception("attr %s not in input df" % col)
        else:

            return self.df[col].to_dict()

    @property
    def sid(self):
        return self.get_attr("sample ID")
    @property
    def R1(self):
        return self.get_attr("R1")
    @property
    def R2(self):
        return self.get_attr("R2")

    def is_somatic(self):
        if len(self.df["Somatic"].unique()) > 1:
            return True
        else:
            return False

    def germline_pair(self):
        t_dict = self.df.to_dict(orient='index')
        for k, d in t_dict.items():
            d["SampleID"] = k
        return t_dict

    def somatic_pair(self):
        """
        :return: nid + tid: {Normal:{info}},{Tumor:{info}}
        """
        pair_dict = {}
        if self.is_somatic():
            for source_name in set(self.df.source_name):
                # todo : use group by
                sub_df = self.df.loc[self.df.source_name == source_name, :]
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
                    pair_dict[key]["source_name"] = source_name
                    pair_dict[key]["Normal"] = normal.loc[nid, :].to_dict()
                    pair_dict[key]["Normal"]["SampleID"] = nid
                    pair_dict[key]["Tumor"] = tumor.loc[tid, :].to_dict()
                    pair_dict[key]["Tumor"]["SampleID"] = tid
            return pair_dict
        else:
            raise Exception("No somatic/pair samples detected")

    def get_full_info(self, odir, gettype='germline'):
        # return formatted/output file name, unify all output here.

        if gettype == 'germline':
            sample_dict = self.germline_pair()
        else:
            sample_dict = self.somatic_pair()

        for sample_id, infodict in sample_dict.items():
            if gettype == 'germline':
                list_infos = [infodict]
            else:
                list_infos = [infodict["Normal"],
                              infodict["Tumor"]]

            for info in list_infos:
                sample_name = info.get("SampleID", '')
                project_name = info.get("project_name", "")
                info["trim_R1"] = os.path.join(config.trim_fmt, "{PE1_id}.clean.fq.gz").format(
                    base=odir,
                    PN=project_name,
                    PE1_id=sample_name + "_R1")
                info["trim_R2"] = os.path.join(config.trim_fmt, "{PE2_id}.clean.fq.gz").format(
                    base=odir,
                    PN=project_name,
                    PE2_id=sample_name + "_R2")
                info["sam"] = config.output_fmt.format(
                    path=odir,
                    PN=project_name,
                    SN=sample_name) + '.sam'
                info["sorted_bam"] = config.output_fmt.format(
                    path=odir,
                    PN=project_name,
                    SN=sample_name) + '_sorted.bam'
                info["recal_bam"] = config.output_fmt.format(
                    path=odir,
                    PN=project_name,
                    SN=sample_name) + '.recal_reads.bam'

                info["sorted_cov"] = config.output_fmt.format(
                    path=odir,
                    PN=project_name,
                    SN=sample_name) + '.sorted_cov.info'
                info["recal_cov"] = config.output_fmt.format(
                    path=odir,
                    PN=project_name,
                    SN=sample_name) + '.recal_cov.info'

        return sample_dict


def validate_df(df, filename):
    template_file = os.path.join(os.path.dirname(__file__),
                                 "data_input.template")
    columns_values = open(template_file).read().strip('\n').split('\t')

    if set(df.columns) != set(columns_values):
        raise Exception("INPUT file has unknown header ."
                        "Should be %s, but %s input" % (";".join(df.columns),
                                                                                     ";".join(columns_values)))

    if df["sample ID"].duplicated().any():
        raise Exception("sample_name has duplicated.")

    chdir = os.path.dirname(os.path.abspath(filename))
    # os.chdir(chdir)
    # print('chdir',chdir)
    for idx, row in df.iterrows():
        # auto implement filepath
        # so easy~~~
        row["R1"] = row["R1"] if pd.isna(row["R1"]) else os.path.join(chdir, row["R1"])
        row["R2"] = row["R2"] if pd.isna(row["R2"]) else os.path.join(chdir, row["R2"])
    return columns_values, df
