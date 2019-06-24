

def modify_vcf():
    cmdline = "vt decompose -s {input_vcf} | vt normalize -r {REF_path} - > {output_vcf}"