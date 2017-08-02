import csv
import os

def get_vcf_info(infile):
    info = {}
    with open(infile) as infh:
        reader = csv.DictReader(infh)
        for rec in reader:
            var_id = (rec['Chr'], rec['Start'], rec['End'], rec['Ref'], rec['Alt'])
            if var_id in info:
                raise Exception('Duplicated var id: %s'%(str(var_id)))
            otherinfo = rec['Otherinfo']
            its = otherinfo.split('\t')
            filter = its[6]
            format = its[-2]
            value = its[-1]
            assert len(format.split(':')) == len(value.split(':'))
            details = dict(zip(format.split(':'), value.split(':')))
            info[var_id] = {'filter': filter,
                            'GT': details.get('GT',''),
                            'AD': details.get('AD',''),
                            'AF': details.get('AF','')}
    return info

samples = ['NY-13_mt2_combined',
 'NYN-7_mt2_low_F',
 'NY-11_mt2_combined',
 'NY-15_mt2_combined',
 'NYN-12_mt2_NotInNormal',
 'NYN-11_mt2_NotInNormal',
 'NYN-14_mt2_low_F',
 'NY-7_mt2_combined',
 'NYN-13_mt2_NotInNormal',
 'NYN-15_mt2_low_F',
 'NYN-15_mt2_NotInNormal',
 'NYN-10_mt2_low_F',
 'NYN-9_mt2_NotInNormal',
 'NYN-12_mt2_low_F',
 'NY-12_mt2_combined',
 'NYN-13_mt2_low_F',
 'NYN-7_mt2_NotInNormal',
 'NY-10_mt2_combined',
 'NYN-9_mt2_low_F',
 'NYN-14_mt2_NotInNormal',
 'NYN-10_mt2_NotInNormal',
 'NY-9_mt2_combined',
 'NY-14_mt2_combined',
 'NYN-11_mt2_low_F']


result_path = '/home/liaoth/project/170123_NY/NY_variant/somatic'
result_path2 = '/home/liaoth/project/170123_NY/NY_variant/somatic/somatic_paired_output'
inpath = '/home/liaoth/project/170123_NY/NY_variant/somatic/output'
outfile = '/home/liaoth/project/170123_NY/NY_variant/somatic/final_output/somatic.result.final.xls'

outfh = open(outfile, 'w')
head = False
for sample in samples:
    #vcf_info_file = os.path.join(result_path,'%s_somatic'%(sample),'%s_mt2.merged.av.anno.csv.hg19_multianno.csv'%(sample))
    if 'combined' in sample:
        vcf_info_file = os.path.join(result_path,
                                 '%s_fastV.merged.anno.hg19_multianno.csv' % (sample))
    else:
        vcf_info_file = os.path.join(result_path2,
                                 '%s.csv' % (sample))
    info = get_vcf_info(vcf_info_file)
    infile = os.path.join(inpath, '%s.filtered.final.v1.xls'%(sample))
    lines = open(infile).readlines()
    if not head:
        header = lines[0].rstrip('\n')
        header = ['sample'] + header.split('\t') + ['filter','GT','AD','AF']
        print >>outfh, '\t'.join(header)
        head = True
    for line in lines[1:]:
        line = line.rstrip('\n')
        its = line.split('\t')
        idx = tuple(its[:5])
        var_info = info[idx]
        outs = [sample] + line.split('\t') + ["%s"%var_info['filter'], "%s"%var_info['GT'],"%s"%var_info['AD'],"%s"%var_info['AF']]
        print >>outfh, '\t'.join(outs)

outfh.close()