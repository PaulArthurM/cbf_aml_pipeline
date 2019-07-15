# A Snakefile for a pre-processing of paired normal/tumour samples.

import json
import re

# Configuration file
configfile: "config.yaml"

# Load json configuration file
CONFIG_JSON = json.load(open(config["SAMPLES"]))

SAMPLES = CONFIG_JSON['samples']


def get_sample_name(sample):
    return re.match("(.+?)\.bam$", sample).group(1)


def get_id(sample):
    return re.search("-(\w+)\.", sample).group(1)

def get_lane(sample):
    return re.search("\.(\d+)$", sample).group(1)

TARGETS = []

BQSR_BAM = []
MERGE = []


for SAMPLE in SAMPLES:
    for TYPE in SAMPLES[SAMPLE]:
        LANES = SAMPLES[SAMPLE][TYPE]
        MERGE.append("/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_{type}-{id}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam".format(sample=SAMPLE, type=TYPE, id=get_id(LANES[0]), lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))




TARGETS.extend(MERGE)


rule all:
    input: TARGETS


# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule mark_duplicates:
    input:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bam"
    output:
        marked_bam="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_marked_duplicates.bam",
        metrics_txt="/data1/scratch/pamesl/projet_cbf/data/metrics/{sample}_marked_dup_metrics.txt"
    shell:
        "gatk MarkDuplicates \
            -I {input} \
            -O {output.marked_bam} \
            -M {output.metrics_txt}"


# Generates recalibration table for Base Quality Score Recalibration (BQSR)
#-L {params.intervals_list} \
#--sequence-dictionary  /data1/scratch/pamesl/projet_cbf/data/hg19_data/reference_hg19/ucsc_hg19.dict \
rule base_recalibrator:
    input:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_marked_duplicates.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample}.table"
    params:
        reference=config["REFERENCE"],
        intervals_list=config["intervals_list"]
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {params.reference} \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/dbsnp_138.b37.vcf \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.b37.vcf \
            -O {output}"


# Apply base quality score recalibration
rule apply_BQSR:
    input:
        table="/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample}.table",
        bam="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_marked_duplicates.bam"
    params:
        reference=config["REFERENCE"]
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_marked_duplicates_BQSR.bam"
    shell:
        "gatk ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"


# Merge multiple sorted alignment files, producing a single sorted output file
rule merge_sam_files:
    input:
        lane_1="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_{type}-{id}.{lane_1}_marked_duplicates_BQSR.bam",
        lane_2="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_{type}-{id}.{lane_2}_marked_duplicates_BQSR.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_{type}-{id}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam"
    shell:
        "gatk MergeSamFiles \
            -I {input.lane_1} \
            -I {input.lane_2} \
            -O {output}"


# Rule for create index from BAM file with samtools index
rule samtools_index:
    input:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{merged_samples}_marked_duplicates_BQSR_merge.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{merged_samples}_marked_duplicates_BQSR_merge.bai"
    shell:
        "samtools index -b {input} {output}"
