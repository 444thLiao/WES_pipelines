import vcf


def fix_vcf4pcgr(input_vcf, output_vcf):
    """
    Because the vcf output from mutect2 have some bad FORMAT
    So we need to fix it and make it suitable for pcgr softwares.
    :param input_vcf:
    :param vcf_bed:
    :return:
    """
    need2pop = ["ALT_F1R2",
                "ALT_F2R1",
                "REF_F1R2",
                "REF_F2R1",
                "QSS"]

    vcf_reader = vcf.Reader(filename=input_vcf)

    for i in need2pop:
        vcf_reader.formats.pop(i, 0)
    vcf_writer = vcf.Writer(open(output_vcf, 'w'),
                            vcf_reader)
    for record in vcf_reader:
        record.FORMAT = ':'.join([_
                                  for _ in str(record.FORMAT).split(':')
                                  if _ not in need2pop])
        for s in record.samples:
            _cache = s.data._asdict()
            for i in need2pop:
                _cache.pop(i, 0)
            new_s_data = vcf_reader._parse_sample_format(record.FORMAT)
            s.data = new_s_data(**_cache)
        vcf_writer.write_record(record)
