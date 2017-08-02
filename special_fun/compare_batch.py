import csv
import os
import shutil


def get_filters_vcf_info(otherinfo):
    # return the filter field from VCF info
    its = otherinfo.split()
    return its[6]


def compare_somatic_to_germline(infile_germline, infile_somatic, outfile_somatic_filtered):
    standard_fields = ['Chr', 'Start', 'End', 'Ref', 'Alt']

    ref_germline = set()

    germline_i = 0
    with open(infile_germline) as infh:
        reader = csv.DictReader(infh)
        for rec in reader:
            idx = [rec[it] for it in standard_fields]
            idx = tuple(idx)
            filter_info = get_filters_vcf_info(rec['Otherinfo'])
            if filter_info == 'PASS':
                ref_germline.add(idx)
            germline_i += 1

    somatic_i = 0
    somatic_filtered_i = 0
    with open(infile_somatic) as infh:
        reader = csv.DictReader(infh)
        field_names = reader.fieldnames

        outfh = open(outfile_somatic_filtered, 'w')
        writer = csv.DictWriter(outfh, fieldnames=field_names)
        writer.writeheader()

        for rec in reader:
            somatic_i += 1
            idx = [rec[it] for it in standard_fields]
            idx = tuple(idx)
            if not idx in ref_germline:
                somatic_filtered_i += 1
                writer.writerow(rec)
        outfh.close()


    summary = {'germline variants':germline_i,
               'somatic variants':somatic_i,
               'somatic variants filtered':somatic_filtered_i}
    return summary

if __name__ == '__main__':
    #samples = ['NYS-1-4-_S13','NYS-2_S14','NYS-3_S15']
    samples = ['NYS-A_S1', 'NYS-B_S2', 'NYS-C_S3', 'NYS-D_S4']

    output_path = '/home/liaoth/project/170123_NY/NYS_variant/output'
    input_path_germline = '/home/liaoth/project/170123_NY/NYS_variant/germline'
    input_path_somatic = '/home/liaoth/project/170123_NY/NYS_variant/somatic'
    result_path = '/home/liaoth/project/170123_NY/NYS_variant/output'
    summary_file = '/home/liaoth/project/170123_NY/NYS_variant/output/AA.summary.xls'

    results = {}
    for sample in samples:
        print sample
        infile_germline = os.path.join(input_path_germline,'%s.merged.anno.hg19_multianno.csv'%(sample,))
        infile_somatic = os.path.join(input_path_somatic,'%s_mt2_fastV.merged.anno.hg19_multianno.csv'%(sample,))
        outfile_somatic_filtered = os.path.join(result_path, '%s_mt2_fastV.merged.anno.hg19_multianno.csv'%(sample,))
        shutil.copy(infile_germline, os.path.join(result_path,'%s.merged.anno.hg19_multianno.csv'%(sample,)))
        summary = compare_somatic_to_germline(infile_germline, infile_somatic, outfile_somatic_filtered)
        results[sample] = summary

    outfh = open(summary_file,'w')
    col_names = ['germline variants', 'somatic variants', 'somatic variants filtered']
    header = ['sample'] + col_names
    print >>outfh, '\t'.join(header)
    for sample in samples:
        outs = [sample]
        for it in col_names:
            outs.append(results[sample][it])
        outs = [str(it) for it in outs]
        print >>outfh, '\t'.join(outs)
    outfh.close()
