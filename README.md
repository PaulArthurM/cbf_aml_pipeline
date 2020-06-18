# Yet Another Pipeline for Normal/Tumor Variant Calling from Whole-Exome Sequencing data v0.1

A Snakemake pipeline for the bioinformatics analysis of WES data with normal/tumor pair.

Multiple variant callers are being added :

- Mutect2 (GATK4)

- Strelka2

- Freebayes

- SomaticSniper


# Genomic data

GRCh37 data have been downloaded with following command-line using [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/):

```
$ aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/ ./references/Homo_sapiens/GATK/GRCh37/
```
