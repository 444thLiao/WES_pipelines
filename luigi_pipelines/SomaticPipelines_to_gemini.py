from collections import defaultdict

import luigi

from luigi_pipelines.SomaticPipelines import MuTect2_single, PrintReads, MuTect2_pair
from luigi_pipelines.share_luigi_tasks import Add_cov_infos, gemini_part, vt_part, vep_part
from special_fun import Add_cov_ino_in_vcf as P_vcf


class Add_cov_infos_SO(Add_cov_infos):
    def requires(self):
        return [MuTect2_single(sample_NT=self.sampleID, dry_run=self.dry_run),
                PrintReads(sampleID=self.sampleID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.mt2.bam', '.mt2.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_SO(self.input()[1].path,
                            self.input()[0].path.replace('.bam', '.vcf'),
                            self.output().path,
                            REF_file_path)


class Add_cov_infos_PA(Add_cov_infos):
    def requires(self):
        _IDs = str(self.sampleID).split(',')
        if pfn(_IDs[0], 'mt2_for') == NORMAL_SIG:
            normal_ID = _IDs[0]
            tumor_ID = _IDs[1]
        else:
            tumor_ID = _IDs[0]
            normal_ID = _IDs[1]

        return [MuTect2_pair(sample_IDs=self.sampleID, dry_run=self.dry_run),
                PrintReads(sampleID=normal_ID, dry_run=self.dry_run),
                PrintReads(sampleID=tumor_ID, dry_run=self.dry_run)]

    def output(self):
        return luigi.LocalTarget(self.input()[0].path.replace('.mt2.bam', '.mt2.added_cov.vcf'))

    def run(self):
        P_vcf.Add_in_vcf_PA([self.input()[1].path, self.input()[2].path],
                            self.input()[0].path.replace('.bam', '.vcf'),
                            self.output().path,
                            REF_file_path)


#########15
class vt_part(vt_part):
    def requires(self):
        if self.sampleID in pair_bucket:
            pair_value = pair_bucket[self.sampleID]
            return [Add_cov_infos_PA(sampleID=','.join(pair_value))]

        else:
            return [Add_cov_infos_SO(sampleID=self.sampleID)]


class vep_part(vep_part):
    def requires(self):
        return vt_part(sampleID=self.sampleID, dry_run=self.dry_run)


class gemini_part(gemini_part):
    def requires(self):
        return vep_part(sampleID=self.sampleID,
                        dry_run=self.dry_run)

#
# class workflow(luigi.Task):
#     x = luigi.Parameter()
#
#     def requires(self):
#         samples_IDs = str(self.x).split(',')
#
#         pair_bucket = defaultdict(list)
#         for _x in samples_IDs:
#             pair_bucket[pfn(_x, 'pair_name')].append(_x)
#         global pair_bucket
#         ###{'XK-2': ['XK-2T_S20', 'XK-2W_S17'],'XK-8': ['XK-8T_S21', 'XK-8W_S18']}
#
#         samples_IDs += [_x for _x in pair_bucket.keys()]
#         for i in samples_IDs:
#             yield new_gemini_part(sampleID=i)

# python -m luigi --module SomaticPipelines_for_NY workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
