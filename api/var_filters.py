#################################################################################
#### Rewrite this script but without changing any logic
####
#### tianhua Liao 2019.06.17
#################################################################################
import sys
from os.path import dirname, join
sys.path.insert(0, dirname(dirname(__file__)))
from parse_file_name import fileparser
from toolkit import valid_path
import re
from collections import defaultdict
import pandas as pd
import click
from tqdm import tqdm

########################################################################################################################
### definition, fields and data binding
########################################################################################################################

standard_fields = ['Chr', 'Start', 'End', 'Ref', 'Alt']
gene_fields = ['Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene']
snp_fields = ['snp138']
clinvar_fields = ['CLINSIG', 'CLNDBN', 'CLNACC', 'CLNDSDB', 'CLNDSDBID']
vcf_fields = ['FILTER', 'QUAL', 'FORMAT', 'genotype']

########################################################################################################################
### filters and functions
########################################################################################################################
########################################################################################################################
### dbsnp common
snp138_common_file = '/home/liaoth/data/humandb/snp138Common.name.txt'
lines = open(snp138_common_file, 'r').readlines()
lines = [it.strip() for it in lines]
snp138_common = set(lines)


def get_filters_vcf_info(otherinfo):
    # return the filter field from VCF info
    its = otherinfo.split()
    return its[6]


def filter_any_common_vars(rec):
    #######################################################
    # filter: snp138_common
    #######################################################
    dbsnp = rec['snp138']
    if (not dbsnp.startswith('rs')) and (not dbsnp == '.'):
        raise Exception('UNKNOWN dbsnp: %s' % (dbsnp))
    if dbsnp.startswith('rs'):
        # if len(dbsnp)!=10:
        #    raise Exception('Unexpected RS_ID: %s'%(dbsnp))
        if dbsnp in snp138_common:
            return False
    #######################################################
    # filter: 1000g2014oct_all
    #######################################################
    if rec['1000g2014oct_all'] != '.':
        if float(rec['1000g2014oct_all']) > 0.01:
            return False
    #######################################################
    # filter: ExAC_ALL
    #######################################################
    if rec['ExAC_ALL'] != '.':
        if float(rec['ExAC_ALL']) > 0.01:
            return False
    #######################################################
    # filter: esp6500siv2_all
    #######################################################
    if rec['esp6500siv2_all'] != '.':
        if float(rec['esp6500siv2_all']) > 0.01:
            return False
    return True


gene_func_filters = ['splicing']
exon_func_filters_rank1 = ['frameshift deletion', 'frameshift insertion', 'stopgain', 'stoploss',
                           'nonframeshift deletion', 'nonframeshift insertion',
                           'nonsynonymous SNV']

exon_func_filters_rank2 = ['frameshift deletion', 'frameshift insertion', 'stopgain', 'stoploss']


def filter_func_vars_rank1(rec):
    if rec['ExonicFunc.refGene'] in exon_func_filters_rank1:
        return True
    if rec['Func.refGene'] in gene_func_filters:
        return True
    return False


def filter_func_vars_rank2(rec):
    if rec['ExonicFunc.refGene'] in exon_func_filters_rank2:
        return True
    return False


def load_infile(infile):
    indf = pd.read_csv(infile, sep=',')
    fields = list(indf.columns)
    if not 'snp138' in fields:
        raise Exception('snp138 need to be annotated with')
    if not '1000g2014oct_all' in fields:
        raise Exception('1000g2014oct_all need to be annotated with')
    if not 'ExAC_ALL' in fields:
        raise Exception('ExAC need to be annotated with')
    if not 'esp6500siv2_all' in fields:
        raise Exception('esp6500siv2_all need to be annotated with')

    return indf


def exec_filter(infile, sample_id, output_prefix):
    outfile_imp = output_prefix + '.imp.xls'
    outfile_patho = output_prefix + '.filtered.patho.xls'
    outfile_final = output_prefix + '.filtered.final.xls'
    ########################################################################################################################
    stats = defaultdict(dict)
    stats[sample_id]['total variants'] = 0
    stats[sample_id]['failed QC variants'] = 0
    stats[sample_id]['dbsnp'] = 0
    stats[sample_id]['dbsnp common'] = 0
    stats[sample_id]['1K genome MAF >0.01'] = 0
    stats[sample_id]['ExAC MAF >0.01'] = 0
    stats[sample_id]['esp6500 MAF >0.01'] = 0
    stats[sample_id]['final variants'] = 0
    stats[sample_id]['clinvar pathogenic variants'] = 0
    stats[sample_id]['important variants'] = 0

    ########################################################################################################################
    columns = standard_fields + gene_fields + snp_fields + clinvar_fields
    imp_df = pd.DataFrame(columns=columns, )
    patho_df = pd.DataFrame(columns=columns, )
    final_df = pd.DataFrame(columns=columns, )

    input_df: pd.DataFrame = load_infile(infile)
    stats[sample_id]['total variants'] = input_df.shape[0]
    for idx, row in tqdm(input_df.iterrows(), total=input_df.shape[0]):

        dbsnp = row['snp138']
        ### clinvar pathogenic records
        clin_sigs = row['CLINSIG']

        ### skip vcf filters failed
        vcf_filter = get_filters_vcf_info(row['Otherinfo'])
        if vcf_filter != 'PASS':
            stats[sample_id]['failed QC variants'] += 1
            continue

        ### report all pathogenic vars
        if 'PATHOGENIC' in clin_sigs.upper():
            outs = []
            outs += [row[k] for k in standard_fields]
            outs += [row[k] for k in gene_fields]
            outs += [row[k] for k in snp_fields]
            clinvar_info = [row[k] for k in clinvar_fields]
            clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
            outs += clinvar_info
            patho_df = patho_df.append(pd.DataFrame([outs], columns=patho_df.columns, index=[0]))

        if (not dbsnp.startswith('rs')) and (not dbsnp == '.'):
            raise Exception('UNKNOWN dbsnp: %s' % (dbsnp))
        if dbsnp.startswith('rs'):
            stats[sample_id]['dbsnp'] += 1
            if dbsnp in snp138_common:
                stats[sample_id]['dbsnp common'] += 1
        if row['1000g2014oct_all'] != '.':
            if float(row['1000g2014oct_all']) > 0.01:
                stats[sample_id]['1K genome MAF >0.01'] += 1
        if row['ExAC_ALL'] != '.':
            if float(row['ExAC_ALL']) > 0.01:
                stats[sample_id]['ExAC MAF >0.01'] += 1
        if row['esp6500siv2_all'] != '.':
            if float(row['esp6500siv2_all']) > 0.01:
                stats[sample_id]['esp6500 MAF >0.01'] += 1

        ### report final vars (rank1)
        if filter_any_common_vars(row) and filter_func_vars_rank1(row):
            outs = []
            outs += [row[k] for k in standard_fields]
            outs += [row[k] for k in gene_fields]
            outs += [row[k] for k in snp_fields]
            clinvar_info = [row[k] for k in clinvar_fields]
            clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
            outs += clinvar_info
            final_df = final_df.append(pd.DataFrame([outs], columns=final_df.columns, index=[0]))

        ### report final vars (rank2)
        if filter_any_common_vars(row) and filter_func_vars_rank2(row):
            outs = []
            outs += [row[k] for k in standard_fields]
            outs += [row[k] for k in gene_fields]
            outs += [row[k] for k in snp_fields]
            clinvar_info = [row[k] for k in clinvar_fields]
            clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
            outs += clinvar_info
            imp_df = imp_df.append(pd.DataFrame([outs], columns=imp_df.columns, index=[0]))

    imp_df.to_csv(outfile_imp, index=False, sep='\t')
    final_df.to_csv(outfile_final, index=False, sep='\t')
    patho_df.to_csv(outfile_patho, index=False, sep='\t')
    stats[sample_id]['final variants'] = final_df.shape[0]
    stats[sample_id]['clinvar pathogenic variants'] = patho_df.shape[0]
    stats[sample_id]['important variants'] = imp_df.shape[0]

    return stats


@click.command()
@click.option("--tab")
@click.option("-i", "--indir")
@click.option("-o", "--outdir")
def main(tab, indir, outdir):
    valid_path(outdir, check_odir=1)
    df = fileparser(tab)
    total_stats = dict()
    project_name = ''
    for idx, row in df.df.iterrows():
        sample_id = idx
        project_name = row["project_name"]
        infile = join(indir, sample_id + '.merged.anno.hg19_multianno.csv')
        stats = exec_filter(infile, sample_id, output_prefix=join(outdir, sample_id))
        total_stats.update(stats)

    summary_df = pd.DataFrame.from_dict(total_stats, orient="index")
    summary_df = summary_df.T
    summary_df.index.name = 'variants'
    summary_df.to_csv(join(outdir, project_name + '.summary.csv'), sep=',', index=1)


if __name__ == '__main__':
    main()
