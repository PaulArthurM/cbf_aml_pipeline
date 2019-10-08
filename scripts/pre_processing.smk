# A Snakefile for a pre-processing of paired normal/tumour samples.


import json
import re
import os.path
from os import listdir
from os.path import isfile, join


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
        if TYPE == 'G':
            LANES = SAMPLES[SAMPLE][TYPE]
            file_1 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[0]), lane=get_lane(LANES[0]))
            file_2 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[1]), lane=get_lane(LANES[1]))
            if (len(LANES)==2):
                if 1:#(os.path.isfile(file_1)) and (os.path.isfile(file_2)):
                    MERGE.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                    MERGE.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            elif (len(LANES)==3):
                file_3 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[2]), lane=get_lane(LANES[2]))
                if 1: #(os.path.isfile(file_1)) and (os.path.isfile(file_2)) and (os.path.isfile(file_3)):
                    MERGE.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    MERGE.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))


#print(MERGE)
TARGETS.extend(MERGE)

rule all:
    input: TARGETS


# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule mark_duplicates:
    input:
        config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}.bam"
    output:
        marked_bam = temp(config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}_marked_duplicates.bam"),
        metrics_txt = config["PROJECT_DIR"] + "data/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
    conda:
        "../envs/gatk4.yaml"
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
        config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
    output:
        temp(config["PROJECT_DIR"] + "data/bam/recal_data_{sample}.table")
    params:
        reference=config["REFERENCE"],
        intervals_list=config["intervals_list"]
    conda:
        "../envs/gatk4.yaml"
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
        table = config["PROJECT_DIR"] + "data/bam/recal_data_{sample}.table",
        bam = config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
    params:
        reference=config["REFERENCE"]
    output:
        temp(config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates_BQSR.bam")
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"


# Merge multiple sorted alignment files, producing a single sorted output file
rule merge_sam_two_files:
    input:
        lane_1 = config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_1}_marked_duplicates_BQSR.bam",
        lane_2 = config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_2}_marked_duplicates_BQSR.bam"
    output:
        config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam"
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk MergeSamFiles \
            -I {input.lane_1} \
            -I {input.lane_2} \
            -O {output}"


# Merge multiple sorted alignment files, producing a single sorted output file
rule merge_sam_three_files:
    input:
        lane_1 = config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_1}_marked_duplicates_BQSR.bam",
        lane_2 = config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_2}_marked_duplicates_BQSR.bam",
        lane_3 = config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_3}_marked_duplicates_BQSR.bam"
    output:
        config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bam"
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk MergeSamFiles \
            -I {input.lane_1} \
            -I {input.lane_2} \
            -I {input.lane_3} \
            -O {output}"


# Rule for create index from BAM file with samtools index
rule samtools_index:
    input:
        config["PROJECT_DIR"] + "data/bam/{merged_samples}_marked_duplicates_BQSR_merge.bam"
    output:
        config["PROJECT_DIR"] + "data/bam/{merged_samples}_marked_duplicates_BQSR_merge.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input} {output}"
