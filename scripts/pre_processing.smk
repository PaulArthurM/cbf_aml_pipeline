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

MERGE_BAM = []
MERGE_BAI = []
FASTQC = []
VCF = []
VCF_IDX = []

for SAMPLE in SAMPLES:
    for TYPE in SAMPLES[SAMPLE]:
        if TYPE == 'G':
            LANES = SAMPLES[SAMPLE][TYPE]
            file_1 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[0]), lane=get_lane(LANES[0]))
            file_2 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[1]), lane=get_lane(LANES[1]))
            if (len(LANES)==2):
                MERGE_BAM.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                MERGE_BAI.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                #FASTQC.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_fastqc.html".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                VCF.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
            elif (len(LANES)==3):
                file_3 = "{project_dir}data/bam/{sample}_{type}-{id}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, id=get_id(LANES[2]), lane=get_lane(LANES[2]))
                MERGE_BAM.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                MERGE_BAI.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                #FASTQC.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_fastqc.html".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                VCF.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
#print(MERGE)
TARGETS.extend(MERGE_BAM)
TARGETS.extend(MERGE_BAI)
#TARGETS.extend(FASTQC)
TARGETS.extend([config["PON_VCF"]])
#TARGETS.extend(VCF)

rule all:
    input: TARGETS


# # Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
# rule mark_duplicates_spark:
#     input:
#         config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}.bam"
#     output:
#         marked_bam = temp(config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}_marked_duplicates.bam"),
#         metrics_txt = config["PROJECT_DIR"] + "data/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
#     conda:
#         "../envs/gatk4.yaml"
#     params:
#         name="mark_duplicates_spark_{sample}_{type}.{lane}",
#         nthread=config["mark_duplicates"]["spark"]["nthread"]
#     shell:
#         "gatk MarkDuplicatesSpark \
#             -I {input} \
#             -O {output.marked_bam} \
#             -M {output.metrics_txt} \
#             --conf 'spark.executor.cores={params.nthread}'"


# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule mark_duplicates:
    input:
        config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}.bam"
    output:
        marked_bam = temp(config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}_marked_duplicates.bam"),
        metrics_txt = config["PROJECT_DIR"] + "data/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
    conda:
        "../envs/gatk4.yaml"
    params:
        name="mark_duplicates_{sample}_{type}.{lane}",
        nthread=config["mark_duplicates"]["classic"]["nthread"]
    shell:
        "gatk MarkDuplicates \
            -I {input} \
            -O {output.marked_bam} \
            -M {output.metrics_txt}"


# rule BQSRPipelineSpark:
#     input:
#         config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
#     output:
#         temp(config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates_BQSR.bam")
#     params:
#         reference=config["reference_GRCh37-lite"],
#         intervals_list=config["intervals_list"],
#         dbsnp_138=config["known-sites"]["dbsnp_138"],
#         mills_1000G=config["known-sites"]["mills_1000G"],
#         name="BQSRPipelineSpark_{sample}",
#         nthread=config["BQSRPipelineSpark"]["nthread"]
#     conda:
#         "../envs/gatk4.yaml"
#     shell:
#         "gatk BQSRPipelineSpark \
#         -R {params.reference} \
#         -I {input} \
#         --known-sites {params.dbsnp_138} \
#         --known-sites {params.mills_1000G} \
#         -L {params.intervals_list} \
#         -O {output} \
#         --conf 'spark.executor.cores={params.nthread}'"


#Generates recalibration table for Base Quality Score Recalibration (BQSR)
rule base_recalibrator:
    input:
        config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
    output:
        temp(config["PROJECT_DIR"] + "data/bam/recal_data_{sample}.table")
    params:
        reference=config["reference_GRCh37-lite"],
        intervals_list=config["intervals_list"],
        name="base_recalibrator_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {params.reference} \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/dbsnp_138.b37.vcf \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.b37.vcf \
            -O {output}"


#Apply base quality score recalibration
rule apply_BQSR:
    input:
        table = config["PROJECT_DIR"] + "data/bam/recal_data_{sample}.table",
        bam = config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
    params:
        reference=config["reference_GRCh37-lite"],
        name="apply_BQSR_{sample}",
        nthread=5
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
    params:
        name="merge_{sample}_{type}.{lane_1}.{lane_2}",
        nthread=config["MergeSamFiles"]["nthread"]
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
    params:
        name="merge_{sample}_{type}.{lane_1}.{lane_2}.{lane_3}",
        nthread=config["MergeSamFiles"]["nthread"]
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
    params:
        name="index_{merged_samples}",
        nthread=config["samtools"]["nthread"]
    shell:
        "samtools index -b {input} {output}"


rule fastqc:
    input:
        config["PROJECT_DIR"] + "data/bam/{merged_samples}_marked_duplicates_BQSR_merge.bam"
    output:
        config["PROJECT_DIR"] + "data/bam/{merged_samples}_marked_duplicates_BQSR_merge_fastqc.html"
    params:
        dir=config["FASTQC"]["DIR"],
        name="fastq_{merged_samples}",
        nthread=1
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -o {params.dir}"


rule Mutect2_tumour_only:
    input:
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bam",
        bai=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bai"
    output:
        temp(config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.gz")
    params:
        ref=config["reference_GRCh37-lite"],
        gnomad=config["mutect2"]["gnomad"]["file"],
        intervals=config["intervals_list"],
        name="Mutect2_tumour_only_{sample}_G.{lane}",
        nthread=config["mutect2"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk Mutect2 \
        -R {params.ref} \
        -I {input.bam} \
        -max-mnp-distance 0 \
        -L {params.intervals} \
        -O {output}"


rule GenomicsDB:
    input:
        #VCF,
        VCF_IDX
    output:
        db=directory(config["db_GDBI"]),
        test="genomicsdb.txt"
    params:
        ref=config["reference_GRCh37-lite"],
        inputString = lambda wildcards, input: " -V ".join(input),
        intervals=config["intervals_list"],
        name="GenomicsDB",
        nthread=20
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk GenomicsDBImport \
        -R {params.ref} \
        -L {params.intervals} \
        --genomicsdb-workspace-path {output.db} \
        --merge-input-intervals \
        -V {params.inputString} && touch {output.test}"


rule CreateSomaticPanelOfNormals:
    input:
        #VCF,
        VCF_IDX,
        test="genomicsdb.txt"
    output:
        config["PON_VCF"]
    params:
        ref=config["reference_GRCh37-lite"],
        db=config["db_GDBI"],
        name="create_PON",
        nthread=20
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk CreateSomaticPanelOfNormals \
        -R {params.ref} \
        -V {params.db} \
        -O {output}"



rule unzip_gz:
    input:
        vcf_gz = config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.gz"
    output:
        vcf = config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf"
    params:
        name="gunzip_{sample}_G.{lane}",
        nthread=1
    shell:
        "gunzip {input.vcf_gz}"



rule IndexFeatureFile:
    input:
        vcf = config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf"
    output:
        vcf_idx = config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.idx"
    params:
        name="IndexFeatureFile_{sample}_G.{lane}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk IndexFeatureFile -F {input.vcf}"
