from api import Annovar2,Annovar1
import luigi







class input_file(luigi.ExternalTask):
    vcf_file = luigi.Parameter()
    def requires(self):
        return luigi.LocalTarget(self.vcf_file)



class new_Annovar1(Annovar1):
    def requires(self):
        return CombineVariants(infodict=self.infodict,
                               dry_run=self.dry_run)


class new_Annovar2(Annovar2):
    def requires(self):
        return [new_Annovar1(infodict=self.infodict,
                            dry_run=self.dry_run)]
