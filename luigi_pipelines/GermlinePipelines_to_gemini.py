###############################################################################################
### Target (Amplicon) sequencing of human exome, germline sample Output a gemini db.
### @GPZ-bioinfo, 20190525
###############################################################################################

import luigi

from luigi_pipelines.GermlinePipelines import CombineVariants, PrintReads
from luigi_pipelines.share_luigi_tasks import luigi_vcf2bed, luigi_bed2info, Add_cov_infos, gemini_part, vt_part, vep_part


class luigi_vcf2bed(luigi_vcf2bed):
    def requires(self):
        return CombineVariants(infodict=self.infodict, dry_run=self.dry_run)


class luigi_bed2info(luigi_bed2info):
    def requires(self):
        return [PrintReads(infodict=self.infodict, dry_run=self.dry_run),
                luigi_vcf2bed(infodict=self.infodict, dry_run=self.dry_run)]


#########14
class Add_cov_infos(Add_cov_infos):

    def requires(self):
        if self.infodict_T is not None:
            tasks = {}
            tasks["Normal"] = [CombineVariants(infodict=self.infodict_N, dry_run=self.dry_run),
                               luigi_bed2info(infodict=self.infodict_N, dry_run=self.dry_run)]
            tasks["Normal"] = [CombineVariants(infodict=self.infodict_T, dry_run=self.dry_run),
                               luigi_bed2info(infodict=self.infodict_T, dry_run=self.dry_run)]
        else:
            return [CombineVariants(infodict=self.infodict_N, dry_run=self.dry_run),
                    luigi_bed2info(infodict=self.infodict_N, dry_run=self.dry_run)]


#########15
class vt_part(vt_part):
    def requires(self):
        return Add_cov_infos(infodict_N=self.infodict, dry_run=self.dry_run)


class vep_part(vep_part):
    def requires(self):
        return vt_part(infodict=self.infodict, dry_run=self.dry_run)


class gemini_part(gemini_part):
    def requires(self):
        return vep_part(infodict=self.infodict,
                        dry_run=self.dry_run)


if __name__ == '__main__':
    luigi.run()

# python -m luigi --module SomaticPipelines_fast_version workflow --x XK-8T_S21,XK-2T_S20,XK-2W_S17,XK-8W_S18 --parallel-scheduling --workers 12 --local-scheduler
