# A Snakefile for a WES analysis pipeline.

import json
import re

# Configuration file
configfile: "config.yaml"

# Load json configuration file
CONFIG_JSON = json.load(open(config["SAMPLES"]))
NORMALS_SAMPLES = json.load(open(config["SAMPLES"]))["panelOfNormals"]

SAMPLES = CONFIG_JSON['samples']
VCF_MAP = config["VCF_MAP"]


def get_sample_name(sample):
    return re.match("(.+?)\.bam$", sample).group(1)


def return_vcf_name(pairs):
    vcf="/data1/scratch/pamesl/projet_cbf/data/vcf/{normal_1}_and{normal_2}_vs_{tumour_1}_and_{tumour_2}_mutect2.vcf"
    return vcf.format(normal_1=pairs[2], normal_2=pairs[3], tumour_1=pairs[0], tumour_2=pairs[1])


#ALL_BAI = []
BQSR_BAM = []
NORMALS = []
VARIANT_CALLING = []
TARGETS = []

for SAMPLE in SAMPLES:
    tumour_1 = SAMPLES[SAMPLE]["D"][0]
    tumour_2 = SAMPLES[SAMPLE]["D"][1]
    normal_1 = SAMPLES[SAMPLE]["G"][0]
    normal_2 = SAMPLES[SAMPLE]["G"][1]
    pairs = [tumour_1, tumour_2, normal_1, normal_2]
    VARIANT_CALLING.append(return_vcf_name(pairs))
    for TYPE in SAMPLES[SAMPLE]:
        for LANE in SAMPLES[SAMPLE][TYPE]:
            sample_name = LANE #get_sample_name(LANE)
            bai_file = "/data1/scratch/pamesl/projet_cbf/data/bam/{sample_name}.bai"
            #ALL_BAI.append(bai_file.format(sample_name=sample_name)) #  Useless for the moment
            bqsr_file = "/data1/scratch/pamesl/projet_cbf/data/bam/{sample_name}_BQSR.bam"
            bqsr_table = "/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample_name}.table"
            BQSR_BAM.append(bqsr_table.format(sample_name=sample_name))
            BQSR_BAM.append(bqsr_file.format(sample_name=sample_name))
    #NORMALS.append()

#for SAMPLE in NORMALS_SAMPLES:
    #TARGETS.append("/data1/scratch/pamesl/projet_cbf/data/vcf/{sample}_single_sample.vcf.gz".format(sample=SAMPLE))


#TARGETS.extend(ALL_BAI)
TARGETS.extend(BQSR_BAM)
TARGETS.extend(VARIANT_CALLING)

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
    params:
        reference=config["REFERENCE"]
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {reference} \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/All_20180423.vcf.gz \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
            -O {output}"


# Apply base quality score recalibration
rule apply_BQSR:
    input:
        table="/data1/scratch/pamesl/projet_cbf/data/bam/recal_data_{sample}.table",
        bam="/data1/scratch/pamesl/projet_cbf/data/bam/{sample}.bam"
    params:
        reference=config["REFERENCE"]
    output:
        "/data1/scratch/pamesl/projet_cbf/data/bam/{sample}_BQSR.bam"
    shell:
        "gatk ApplyBQSR \
            -R {reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"


rule variant_calling_Mutect2:
    input:
        tumour_bam_1="/data1/scratch/pamesl/projet_cbf/data/bam/{tumour_1}_BQSR.bam",
        tumour_bam_2="/data1/scratch/pamesl/projet_cbf/data/bam/{tumour_2}_BQSR.bam",
        normal_bam_1="/data1/scratch/pamesl/projet_cbf/data/bam/{normal_1}_BQSR.bam",
        normal_bam_2="/data1/scratch/pamesl/projet_cbf/data/bam/{normal_2}_BQSR.bam",
        pon=config['PON_VCF']
    output:
        "/data1/scratch/pamesl/projet_cbf/data/vcf/{normal_1}_and{normal_2}_vs_{tumour_1}_and_{tumour_2}_mutect2.vcf"
    params:
        reference=config["REFERENCE"]
    shell:
        "gatk Mutect2 \
            -R {reference} \
            -I {input.tumour_bam_1}} \
            -I {input.tumour_bam_2} \
            -I {input.normal_bam_1} \
            -I {input.normal_bam_2} \
            -normal {normal_1} \
            -normal {normal_2} \
            --germline-resource af-only-gnomad.vcf.gz \
            --panel-of-normals {input.pon} \
            -O {normal_1}_{normal_2}_vs_{tumour_1}_and_{tumour_2}.vcf.gz"


rule create_vcf_for_normal:
    input:
        normal="/data1/scratch/pamesl/projet_cbf/data/bam/{normal}.bam"
    params:
        reference=config["REFERENCE"]
    output:
        vcf="/data1/scratch/pamesl/projet_cbf/data/vcf/{normal}_single_sample.vcf.gz",
        done=touch("create_vcf_for_normal.done")
    shell:
        "gatk Mutect2 \
            -R {reference} \
            -I {input.normal} \
            -O {output.vcf}"


rule create_DB_GenomicsDBImport:
    input:
        #test="create_vcf_for_normal.done"
        expand("/data1/scratch/pamesl/projet_cbf/data/vcf/{normal}_single_sample.vcf.gz", normal=NORMALS_SAMPLES)
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
        #db=directory(config["db_GDBI"]),
        reference=config["REFERENCE"]
    shell:
        "gatk CreateSomaticPanelOfNormals \
        -R {params.reference} \
        -V gendb://{input.db} \
        -O {output.pon}"
