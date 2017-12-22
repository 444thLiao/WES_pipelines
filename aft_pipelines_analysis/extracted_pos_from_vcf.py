import os


def prepare_vcf(vcf_path):
    if not vcf_path.endswith('.gz'):
        os.system('`which bgzip` %s' % vcf_path)
        vcf_path += '.gz'
        print vcf_path
        os.system('`which bcftools` index %s' % vcf_path)
    else:
        if not os.path.isfile(vcf_path+'.csi'):
            os.system('`which bcftools` index %s' % vcf_path)

def merge_two_vcf(pair_vcf,bed,single_vcf,output_vcf):
    prepare_vcf(pair_vcf)
    prepare_vcf(single_vcf)
    pair_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % pair_vcf).read().replace('\n','').replace('\t',',')
    single_sample_name = os.popen("zgrep '^#C' %s | cut -f 10-" % single_vcf).read().replace('\n','').replace('\t',',')

    formatted_line = '''bcftools view {vcf} -R {bed} -s ^{SM} --force-samples | vcf-sort -c > {output}'''
    cmdline1 = formatted_line.format(vcf=pair_vcf+'.gz',bed=bed,output=output_vcf+'1',SM=pair_sample_name)
    cmdline2 = formatted_line.format(vcf=single_vcf+'.gz',bed=bed,output=output_vcf+'2',SM=single_sample_name)
    print cmdline1,cmdline2
    os.system(cmdline1)
    os.system(cmdline2)

    prepare_vcf(output_vcf+'1');prepare_vcf(output_vcf+'2');

    formatted_line2 = """bcftools concat {o1} {o2} -a -d all > {output}""".format(o1=output_vcf+'1.gz',o2=output_vcf+'2.gz',output=output_vcf)
    os.system(formatted_line2)
    os.system('rm {o1}* ;rm {o2}*'.format(o1=output_vcf + '1', o2=output_vcf + '2'))
    prepare_vcf(output_vcf)




if __name__ =='__main__':
    pass