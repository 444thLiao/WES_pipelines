import os
cmdline = "python ~/tools/pcgr-0.3.4/pcgr.py --input_vcf '{ivcf}' ~/tools/pcgr-0.3.4/ '{odir}' {oname}"

cmdline2 = "python ~/tools/pcgr-0.5.3/pcgr.py --input_vcf '{ivcf}' ~/tools/pcgr-0.5.3/ '{odir}' {oname}"
if __name__ == '__main__':
    ref = '/home/liaoth/data/hg19/ucsc.hg19.fasta'

    ivcf = '/home/liaoth/project/180104_XK/gpz_server/vcf_storge/XK-57_TS_merged.mt2.vcf.gz'
    odir = '/home/liaoth/project/180104_XK/gpz_server/gemini_annotate/pcgr/XK-32_FINAL'
    oname = 'XK-32_FINAL'

    print cmdline.format(ivcf=ivcf,odir=odir,oname=oname)