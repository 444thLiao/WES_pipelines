import os


def prepare_vcf(vcf_path):
    if not vcf_path.endswith('.gz'):
        os.system('`which bgzip` %s' % vcf_path)
        vcf_path += '.gz'
        print vcf_path
        os.system('/home/liaoth/tools/bcftools-1.3.1/bcftools index %s' % vcf_path)
    else:
        if not os.path.isfile(vcf_path + '.csi'):
            os.system('/home/liaoth/tools/bcftools-1.3.1/bcftools index %s' % vcf_path)


def merge_two_vcf(pair_vcf, bed, single_vcf, output_vcf):
    prepare_vcf(pair_vcf)
    prepare_vcf(single_vcf)
    pair_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % pair_vcf).read().replace('\n', '').replace('\t', ',')
    single_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % single_vcf).read().replace('\n', '').replace('\t',
                                                                                                               ',')

    formatted_line = '''/home/liaoth/tools/bcftools-1.3.1/bcftools view {vcf} -R {bed} -s ^{SM} --force-samples | vcf-sort > {output}'''
    cmdline1 = formatted_line.format(vcf=pair_vcf.replace('.gz', '') + '.gz', bed=bed, output=output_vcf + '1',
                                     SM=pair_sample_name)
    cmdline2 = formatted_line.format(vcf=single_vcf.replace('.gz', '') + '.gz', bed=bed, output=output_vcf + '2',
                                     SM=single_sample_name)
    print cmdline1, cmdline2
    os.system(cmdline1)
    os.system(cmdline2)

    prepare_vcf(output_vcf + '1');
    prepare_vcf(output_vcf + '2');

    formatted_line2 = """/home/liaoth/tools/bcftools-1.3.1/bcftools concat {o1} {o2} -a -d all > {output}""".format(
        o1=output_vcf + '1.gz',
        o2=output_vcf + '2.gz',
        output=output_vcf)
    os.system(formatted_line2)
    os.system('rm {o1}* ;rm {o2}*'.format(o1=output_vcf + '1', o2=output_vcf + '2'))
    prepare_vcf(output_vcf)


if __name__ == '__main__':
    import sys

    pairedvcf = sys.argv[1]
    bed_f = sys.argv[2]
    single_vcf = sys.argv[3]
    output_vcf = sys.argv[4]
    # merge_two_vcf('/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25.mt2.vcf',
    #               '/home/liaoth/project/180104_XK/gpz_server/analysis_result/filtered_somatic/XK-25_filtered_somatic.bed'
    #               , '/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25T.mt2.vcf',
    #               '/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-25_TS_merged.mt2.vcf')
    merge_two_vcf(pairedvcf, bed_f, single_vcf, output_vcf)
