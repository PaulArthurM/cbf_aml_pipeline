# A Snakefile for a WES analysis pipeline.

import json
import re

# Configuration file
configfile: "config.yaml"

# Load json configuration file
CONFIG_JSON = json.load(open(config["SAMPLES"]))

SAMPLES = CONFIG_JSON['samples']


def get_sample_name(sample):
    return re.match("(.+?)\.bam$", sample).group(1)


ALL_BAI = []
BQSR_BAM = []
NORMALS = []
TARGETS = []

for SAMPLE in SAMPLES:
    for TYPE in SAMPLES[SAMPLE]:
        print(SAMPLES[SAMPLE][TYPE])
        for LANE in SAMPLES[SAMPLE][TYPE]:
            print(LANE)
            sample_name = LANE #get_sample_name(LANE)
            bai_file = "/data1/scratch/pamesl/projet_cbf/data/bam/{sample_name}.bai"
            ALL_BAI.append(bai_file.format(sample_name=sample_name))
            bqsr_file = "/data1/scratch/pamesl/projet_cbf/data/bam/{sample_name}_BQSR.bam"
            bqsr_table = "/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample_name}.table"
            BQSR_BAM.append(bqsr_table.format(sample_name=sample_name))
            BQSR_BAM.append(bqsr_file.format(sample_name=sample_name))
    #NORMALS.append()


TARGETS.extend(ALL_BAI)
TARGETS.extend(BQSR_BAM)

rule all:
    input: TARGETS


# Rule for create index from BAM file with samtools index
rule samtools_index:
    input:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bai"
    shell:
        "samtools index -b {input} {output}"


# Merge multiple sorted alignment files, producing a single sorted output file
#rule merge_sam_files:
#    input:



# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule mark_duplicates:
    input:
        "data/bam/{sample}.bam"
    output:
        marked_bam="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_marked_duplicates.bam",
        metrics_txt="/data1/scratch/pamesl/projet_cbf/data/metrics/{sample}_marked_dup_metrics.txt"
    shell:
        "conda activate gatk4_4.1.2.0_env &&"
        "java -jar picard.jar MarkDuplicates \
            I={input} \
            O={output.marked_bam} \
            M={output.metrics_txt} &&"
        "conda deactivate"


# Generates recalibration table for Base Quality Score Recalibration (BQSR)
rule base_recalibrator:
    input:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample}.table"
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R reference.fasta \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/All_20180423.vcf.gz \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
            -O {output}"


# Apply base quality score recalibration
rule apply_BQSR:
    input:
        table="/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample}.table",
        bam="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_BQSR.bam"
    shell:
        "gatk ApplyBQSR \
            -R reference.fasta \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"


rule variant_calling_Mutect2:
    input:
        tumour_bam="-I /data1/scratch/pamesl/projet_cbf/data/bam/{tumour_1}.bam -I /data1/scratch/pamesl/projet_cbf/data/bam/{tumour_2}.bam",
        normal_bam="-I /data1/scratch/pamesl/projet_cbf/data/bam/{normal_1}.bam -I /data1/scratch/pamesl/projet_cbf/data/bam/{normal_2}.bam"
    output:
        "/data1/scratch/pamesl/projet_cbf/data/vcf/{normal_1}_and{normal_2}_vs_{tumour_1}_and_{tumour_2}_mutect2.vcf"
    shell:
        "gatk Mutect2 \
            -R reference.fa \
            {input.tumour_bam}} \
            {input.normal_bam} \
            -normal {normal_1} \
            -normal {normal_2} \
            --germline-resource af-only-gnomad.vcf.gz \
            --panel-of-normals pon.vcf.gz \
            -O {normal_1}_{normal_2}_vs_{tumour_1}_and_{tumour_2}.vcf.gz"
