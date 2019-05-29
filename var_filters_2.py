import csv
import glob
import re
from collections import defaultdict

# from av_info import get_av_info

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


########################################################################################################################
def Combined_analysis(args, input_path, output_path, input_path_format):
    ########################################################################################################################
    ### dbsnp common
    lines = open(snp138_common_file, 'r').readlines()
    lines = [it.strip() for it in lines]
    snp138_common = set(lines)

    ### samples
    samples = args.split(',')
    results = defaultdict(dict)

    outfile_summary = "%s/vars_summary.v1.xls" % output_path

    for sample in samples:
        print sample
        infile = '%s/%s' % (input_path, input_path_format.format(SN=sample))
        outfile_imp = "%s/%s.imp.v1.xls" % (output_path, sample,)
        outfile_patho = "%s/%s.filtered.patho.v1.xls" % (output_path, sample,)
        outfile_final = "%s/%s.filtered.final.v1.xls" % (output_path, sample,)

        ########################################################################################################################
        total_vars = 0
        dbsnp_vars = 0
        dbsnp_common_vars = 0
        one_thousand_maf = 0
        exac_maf = 0
        esp6500_maf = 0
        patho_vars = 0
        important_vars = 0
        patho_vars = 0
        failed_qc_vars = 0
        final_vars = 0
        ########################################################################################################################

        i = 0

        outfh_imp = open(outfile_imp, 'w')
        outfh_patho = open(outfile_patho, 'w')
        outfh_final = open(outfile_final, 'w')

        # out_fields = standard_fields + gene_fields + vcf_fields + snp_fields + clinvar_fields
        out_fields = standard_fields + gene_fields + snp_fields + clinvar_fields
        print >> outfh_imp, '\t'.join(out_fields)
        print >> outfh_patho, '\t'.join(out_fields)
        print >> outfh_final, '\t'.join(out_fields)

        with open(infile) as infh:
            reader = csv.DictReader(infh)
            fields = reader.fieldnames
            ### check for data binding
            if not 'snp138' in fields:
                raise Exception('snp138 need to be annotated with')
            if not '1000g2014oct_all' in fields:
                raise Exception('1000g2014oct_all need to be annotated with')
            if not 'ExAC_ALL' in fields:
                raise Exception('ExAC need to be annotated with')
            if not 'esp6500siv2_all' in fields:
                raise Exception('esp6500siv2_all need to be annotated with')

            for rec in reader:
                total_vars += 1
                var_id = (rec['Chr'],
                          rec['Start'],
                          rec['End'],
                          rec['Ref'],
                          rec['Alt'])
                dbsnp = rec['snp138']
                ### clinvar pathogenic records
                clin_sigs = rec['CLINSIG']

                ### skip vcf filters failed
                vcf_filter = get_filters_vcf_info(rec['Otherinfo'])
                if vcf_filter != 'PASS':
                    failed_qc_vars += 1
                    continue

                ### report all pathogenic vars
                if 'PATHOGENIC' in clin_sigs.upper():
                    patho_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_patho, '\t'.join(outs)

                if (not dbsnp.startswith('rs')) and (not dbsnp == '.'):
                    raise Exception('UNKNOWN dbsnp: %s' % (dbsnp))
                if dbsnp.startswith('rs'):
                    # if len(dbsnp)!=10:
                    #    raise Exception('Unexpected RS_ID: %s'%(dbsnp))
                    dbsnp_vars += 1
                    if dbsnp in snp138_common:
                        dbsnp_common_vars += 1
                if rec['1000g2014oct_all'] != '.':
                    if float(rec['1000g2014oct_all']) > 0.01:
                        one_thousand_maf += 1
                if rec['ExAC_ALL'] != '.':
                    if float(rec['ExAC_ALL']) > 0.01:
                        exac_maf += 1
                if rec['esp6500siv2_all'] != '.':
                    if float(rec['esp6500siv2_all']) > 0.01:
                        esp6500_maf += 1

                ### report final vars (rank1)
                if filter_any_common_vars(rec) and filter_func_vars_rank1(rec):
                    final_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_final, '\t'.join(outs)

                ### report final vars (rank2)
                if filter_any_common_vars(rec) and filter_func_vars_rank2(rec):
                    important_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_imp, '\t'.join(outs)

        outfh_patho.close()
        outfh_final.close()
        outfh_imp.close()

        results['total variants'][sample] = total_vars
        results['failed QC variants'][sample] = failed_qc_vars
        results['dbsnp'][sample] = dbsnp_vars
        results['dbsnp common'][sample] = dbsnp_common_vars
        results['1K genome MAF >0.01'][sample] = one_thousand_maf
        results['ExAC MAF >0.01'][sample] = exac_maf
        results['esp6500 MAF >0.01'][sample] = esp6500_maf
        results['final variants'][sample] = final_vars
        results['clinvar pathogenic variants'][sample] = patho_vars
        results['important variants'][sample] = important_vars

    ### summary of all samples

    outfh_summary = open(outfile_summary, 'w')
    header = ['variants'] + samples
    print >> outfh_summary, '\t'.join(header)
    row_names = ['total variants', 'failed QC variants',
                 'dbsnp', 'dbsnp common',
                 '1K genome MAF >0.01', 'ExAC MAF >0.01', 'esp6500 MAF >0.01',
                 'important variants', 'clinvar pathogenic variants', 'final variants']
    for row_name in row_names:
        outs = [row_name]
        for sample in samples:
            outs.append(str(results[row_name][sample]))
        print >> outfh_summary, '\t'.join(outs)
    outfh_summary.close()


if __name__ == '__main__':

    ########################################################################################################################
    ### dbsnp common
    snp138_common_file = '/home/liaoth/data/humandb/snp138Common.name.txt'
    lines = open(snp138_common_file, 'r').readlines()
    lines = [it.strip() for it in lines]
    snp138_common = set(lines)

    ### samples
    samples = ['HCC591L', 'HCC591T']
    results = defaultdict(dict)

    outfile_summary = '/home/liaoth/project/YZ_XR_WES/server_result/somatic/YZ_XR_WES_germline_result.xls'

    for sample in samples:
        print sample
        infile = glob.glob(
            '/home/liaoth/project/YZ_XR_WES/server_result/germline/%s.merged.anno.hg19_multianno.csv' % sample)[0]
        outfile_imp = infile.partition('.')[0] + '.imp.v1.xls'
        outfile_patho = infile.partition('.')[0] + '.filtered.patho.v1.xls'
        outfile_final = infile.partition('.')[0] + '.filtered.final.v1.xls'

        ########################################################################################################################
        total_vars = 0
        dbsnp_vars = 0
        dbsnp_common_vars = 0
        one_thousand_maf = 0
        exac_maf = 0
        esp6500_maf = 0
        patho_vars = 0
        important_vars = 0
        patho_vars = 0
        failed_qc_vars = 0
        final_vars = 0
        ########################################################################################################################

        i = 0

        outfh_imp = open(outfile_imp, 'w')
        outfh_patho = open(outfile_patho, 'w')
        outfh_final = open(outfile_final, 'w')

        # out_fields = standard_fields + gene_fields + vcf_fields + snp_fields + clinvar_fields
        out_fields = standard_fields + gene_fields + snp_fields + clinvar_fields
        print >> outfh_imp, '\t'.join(out_fields)
        print >> outfh_patho, '\t'.join(out_fields)
        print >> outfh_final, '\t'.join(out_fields)

        with open(infile) as infh:
            reader = csv.DictReader(infh)
            fields = reader.fieldnames
            ### check for data binding
            if not 'snp138' in fields:
                raise Exception('snp138 need to be annotated with')
            if not '1000g2014oct_all' in fields:
                raise Exception('1000g2014oct_all need to be annotated with')
            if not 'ExAC_ALL' in fields:
                raise Exception('ExAC need to be annotated with')
            if not 'esp6500siv2_all' in fields:
                raise Exception('esp6500siv2_all need to be annotated with')

            for rec in reader:
                total_vars += 1
                var_id = (rec['Chr'], rec['Start'], rec['End'], rec['Ref'], rec['Alt'])
                dbsnp = rec['snp138']
                ### clinvar pathogenic records
                clin_sigs = rec['CLINSIG']

                ### skip vcf filters failed
                vcf_filter = get_filters_vcf_info(rec['Otherinfo'])
                if vcf_filter != 'PASS':
                    failed_qc_vars += 1
                    continue

                ### report all pathogenic vars
                if 'PATHOGENIC' in clin_sigs.upper():
                    patho_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_patho, '\t'.join(outs)

                if (not dbsnp.startswith('rs')) and (not dbsnp == '.'):
                    raise Exception('UNKNOWN dbsnp: %s' % (dbsnp))
                if dbsnp.startswith('rs'):
                    # if len(dbsnp)!=10:
                    #    raise Exception('Unexpected RS_ID: %s'%(dbsnp))
                    dbsnp_vars += 1
                    if dbsnp in snp138_common:
                        dbsnp_common_vars += 1
                if rec['1000g2014oct_all'] != '.':
                    if float(rec['1000g2014oct_all']) > 0.01:
                        one_thousand_maf += 1
                if rec['ExAC_ALL'] != '.':
                    if float(rec['ExAC_ALL']) > 0.01:
                        exac_maf += 1
                if rec['esp6500siv2_all'] != '.':
                    if float(rec['esp6500siv2_all']) > 0.01:
                        esp6500_maf += 1

                ### report final vars (rank1)
                if filter_any_common_vars(rec) and filter_func_vars_rank1(rec):
                    final_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_final, '\t'.join(outs)

                ### report final vars (rank2)
                if filter_any_common_vars(rec) and filter_func_vars_rank2(rec):
                    important_vars += 1
                    outs = []
                    outs += [rec[k] for k in standard_fields]
                    outs += [rec[k] for k in gene_fields]
                    outs += [rec[k] for k in snp_fields]
                    clinvar_info = [rec[k] for k in clinvar_fields]
                    clinvar_info = [re.sub(r'\\x2c', ',', it) for it in clinvar_info]
                    outs += clinvar_info
                    print >> outfh_imp, '\t'.join(outs)

        outfh_patho.close()
        outfh_final.close()
        outfh_imp.close()

        results['total variants'][sample] = total_vars
        results['failed QC variants'][sample] = failed_qc_vars
        results['dbsnp'][sample] = dbsnp_vars
        results['dbsnp common'][sample] = dbsnp_common_vars
        results['1K genome MAF >0.01'][sample] = one_thousand_maf
        results['ExAC MAF >0.01'][sample] = exac_maf
        results['esp6500 MAF >0.01'][sample] = esp6500_maf
        results['final variants'][sample] = final_vars
        results['clinvar pathogenic variants'][sample] = patho_vars
        results['important variants'][sample] = important_vars

    ### summary of all samples

    outfh_summary = open(outfile_summary, 'w')
    header = ['variants'] + samples
    print >> outfh_summary, '\t'.join(header)
    row_names = ['total variants',
                 'failed QC variants',
                 'dbsnp',
                 'dbsnp common',
                 '1K genome MAF >0.01',
                 'ExAC MAF >0.01',
                 'esp6500 MAF >0.01',
                 'important variants',
                 'clinvar pathogenic variants',
                 'final variants']
    for row_name in row_names:
        outs = [row_name]
        for sample in samples:
            outs.append(str(results[row_name][sample]))
        print >> outfh_summary, '\t'.join(outs)

    outfh_summary.close()
