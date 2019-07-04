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


rule create_somatic_panelOfNormals:
    shell:
        "gatk GenomicsDBImport \
            --genomicsdb-workspace-path /data1/scratch/pamesl/projet_cbf/data/GenomicsDBImport \
            --batch-size 50 \
            -L 500 \
            --sample-name-map cohort.sample_map \
            --reader-threads 5"
