# A Snakefile for build Panel Of Normals


import re
import json


# Configuration file
configfile: "config.yaml"


# Load json configuration file
NORMALS_SAMPLES = json.load(open(config["SAMPLES"]))["panelOfNormals"]


VCF_MAP = config["VCF_MAP"]


TARGETS = []
TARGETS.append(config['PON_VCF'])
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
        db=directory(config["db_GDBI"])
    params:
        sample_map=config["VCF_MAP"],
        reader_threads=config["reader_threads"],
        batch_size=config["batch_size"],
        intervals_list=config["intervals_list"],
        reference=config["REFERENCE"]
    shell:
        "gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.db} \
            --batch-size {params.batch_size} \
            -L {params.intervals_list} \
            --sample-name-map {params.sample_map} \
            --reader-threads {params.reader_threads} \
            -R {params.reference}"


# Wait to see GenomicsDBImport output files
rule create_somatic_panelOfNormals:
    input:
        db=config["db_GDBI"]
    output:
        pon=config["VCF_MAP"]
    params:
        db=directory(config["db_GDBI"]),
        reference=config["REFERENCE"]
    shell:
        "gatk CreateSomaticPanelOfNormals \
        -R {params.reference} \
        -V gendb://{params.db} \
        -O {output.pon}"
