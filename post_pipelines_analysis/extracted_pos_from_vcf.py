import sys
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
import os
from luigi_pipelines import config, run_cmd

bcftools_path = config.bcftools_path


# version : 1.9-183-ga5eebf8
def prepare_vcf(vcf_path, log_file=sys.stdout):
    if not vcf_path.endswith('.gz'):
        run_cmd('{bgzip} -c {vcf_p} > {vcf_p}.gz'.format(bgzip=config.bgzip_pro,
                                                         vcf_p=vcf_path),
                log_file=log_file)
        vcf_path += '.gz'
        run_cmd('%s index %s' % (bcftools_path, vcf_path),
                log_file=log_file)
    else:
        if not os.path.isfile(vcf_path + '.csi'):
            run_cmd('%s index %s' % (bcftools_path, vcf_path),
                    log_file=log_file)
    return vcf_path


def merge_two_vcf(pair_vcf, single_vcf, bed, output_vcf, log_file=sys.stdout):
    prepare_vcf(pair_vcf,log_file=log_file)
    prepare_vcf(single_vcf,log_file=log_file)
    pair_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % pair_vcf).read().replace('\n', '').replace('\t', ',')
    single_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % single_vcf).read().replace('\n', '').replace('\t',
                                                                                                               ',')

    formatted_line = '''{bcftools_path} view {vcf} -R {bed} -s ^{SM} --force-samples | {bcftools_path} sort > {output}; '''
    cmdline1 = formatted_line.format(vcf=pair_vcf.replace('.gz', '') + '.gz',
                                     bcftools_path=bcftools_path,
                                     bed=bed,
                                     output=output_vcf + '1',
                                     SM=pair_sample_name)
    cmdline2 = formatted_line.format(vcf=single_vcf.replace('.gz', '') + '.gz',
                                     bcftools_path=bcftools_path,
                                     bed=bed,
                                     output=output_vcf + '2',
                                     SM=single_sample_name)
    run_cmd(cmdline1, log_file=log_file)
    run_cmd(cmdline2, log_file=log_file)

    prepare_vcf(output_vcf + '1',log_file=log_file)
    prepare_vcf(output_vcf + '2',log_file=log_file)

    formatted_line2 = """{bcftools} concat {o1} {o2} -a -d all > {output}""".format(
        bcftools=bcftools_path,
        o1=output_vcf + '1.gz',
        o2=output_vcf + '2.gz',
        output=output_vcf+'3')
    run_cmd(formatted_line2, log_file=log_file)


    with open(output_vcf+'3',"r") as fr:
        with open(output_vcf,'w') as f1:
            for row in fr:
                if row.startswith("##SAMPLE=<ID="):
                    continue
                f1.write(row)
    prepare_vcf(output_vcf, log_file=log_file)
    run_cmd('rm {o}1* ;rm {o}2* ; rm {o}3* '.format(o=output_vcf))

if __name__ == '__main__':
    import sys

    setting_f = sys.argv[1]
    bed_f = sys.argv[2]
    pairedvcf = sys.argv[3]
    single_vcf = sys.argv[4]
    output_vcf = sys.argv[5]
    # merge_two_vcf('/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25.mt2.vcf',
    #               '/home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-25_filtered_somatic.bed'
    #               , '/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25T.mt2.vcf',
    #               '/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25_TS_merged.mt2.vcf')
    merge_two_vcf(pairedvcf, bed_f, single_vcf, output_vcf)
