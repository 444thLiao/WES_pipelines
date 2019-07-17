# WES pipelines

This pipelines is mainly designed for WES(Whole Exome Sequencing) of human, particularly the paired cancer samples.

Including two pipelines(germline & somatic) 


# installation

```bash
git clone https://github.com/444thLiao/WES_pipelines.git
cd WES_pipelines
pip install -r requirement.txt
```

others requirement need to manually install.
* [bcftools](http://samtools.github.io/bcftools/)
* [vt](https://github.com/atks/vt.git)
* [vep](http://asia.ensembl.org/info/docs/tools/vep/index.html)
* [annovar](http://annovar.openbioinformatics.org/en/latest/)
* [gatk](https://software.broadinstitute.org/gatk/download/)
* [samtools](http://samtools.sourceforge.net/)
* [bgzip](http://www.htslib.org/doc/tabix.html)
* [tabix](http://www.htslib.org/doc/tabix.html)
* [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
* [picard](https://broadinstitute.github.io/picard/)
* [pcgr](https://github.com/sigven/pcgr)

some packages are needed because of `pcgr` and its pipelines. So if you want to output result and directly output pcgr-like output. You need to install `pcgr` `vep` and `vt`.

# main pipelines description

After clone this repository, please first check the `setting.py` file and verify the path is existed at your personal computer.

 1. Download  required WES data
 2. Fulfil a table according to `data_input.template`
 2. perform analysis with `main.py`

With below command,

```
python3 WES_pipelines/luigi_pipelines/main.py run -- main_entry --tab WES_pipelines/test_set/somatic/data_input.tsv --odir ~/temp/test_somatic --analysis-type somatic --workers 5 --log-path ~/temp/test_somatic/cmd_log.txt
```
Above code is performing pipelines used for `somatic variants calling`. Finally, it will output to `~/temp/test_somatic`

If you want to perform `germline variants calling` pipelines, you just need to change the parameter of `analysis-type` into `germline`

Validate parameter of `analysis-type` includes:

* germline
* germline_gatk4
* somatic
* somatic_gatk4


# testing data

If you only want to test the completeness of all required software, you could use below command which use testing data to validate.
```
python3 WES_pipelines/luigi_pipelines/main.py test-germline -o ~/temp/test_germline
```

If you want to test others pipelines, you could also change the parameter of `test-germline` into `test-somatic`