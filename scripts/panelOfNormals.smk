# A Snakefile for build Panel Of Normals


import re
import json


# Configuration file
configfile: "config.yaml"


# Load json configuration file
NORMALS_SAMPLES = json.load(open(config["SAMPLES"]))["panelOfNormals"]


VCF_MAP = config["VCF_MAP"]


TARGETS = []
for SAMPLE in NORMALS_SAMPLES:
    TARGETS.append("/data1/scratch/pamesl/projet_cbf/data/vcf/{sample}_single_sample.vcf.gz".format(sample=SAMPLE))


rule all:
    input: TARGETS


rule create_vcf_for_normal:
    input:
        normal="/data1/scratch/pamesl/projet_cbf/data/bam/{normal}.bam"
    output:
        vcf="/data1/scratch/pamesl/projet_cbf/data/vcf/{normal}_single_sample.vcf.gz"
    shell:
        "gatk Mutect2 \
            -R reference.fa \
            -I {input.normal} \
            -O {output.vcf}"


rule create_DB_GenomicsDBImport:
    input:
        test="/data1/scratch/pamesl/projet_cbf/data/vcf/EGAR00001347180_SJCBF016_G-C0DG1ACXX.5_single_sample.vcf.gz"
    output:
        directory("/data1/scratch/pamesl/projet_cbf/data/GenomicsDBImport")
    params:
        sample_map= config["VCF_MAP"],
        reader_threads=config["reader_threads"],
        batch_size=config["batch_size"],
        intervals_size=config["intervals_size"]
    shell:
        "gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output} \
            --batch-size {params.batch_size} \
            -L {params.intervals_size} \
            --sample-name-map {params.sample_map} \
            --reader-threads {params.reader_threads} \
            -R reference.fa"


rule create_somatic_panelOfNormals:
