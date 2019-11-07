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
VCF_SOMATIC = []
VCF_FILERED = []

for SAMPLE in SAMPLES:
    for TYPE in SAMPLES[SAMPLE]:
        if True: #if TYPE == 'D':
            LANES = SAMPLES[SAMPLE][TYPE]
            file_1 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[0]))
            file_2 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[1]))
            if ( (len(LANES)==2) and (os.path.isfile(file_1)) and (os.path.isfile(file_2)) ):
                print(file_1)
                print(file_2)
                MERGE_BAM.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                MERGE_BAI.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                FASTQC.append("{project_dir}{fastq_dir}{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_fastqc.html".format(fastq_dir=config["FASTQC"]["DIR"], project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                VCF.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                #VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1])))
                vcf_somatic = "{project_dir}data/vcf/{sample}_{lanes_normal}-{lanes_tumour}_somatic.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=".".join([get_lane(SAMPLES[SAMPLE]['G'][0]), get_lane(SAMPLES[SAMPLE]['G'][1])]), lanes_tumour=".".join([get_lane(SAMPLES[SAMPLE]['D'][0]), get_lane(SAMPLES[SAMPLE]['D'][1])]))
                vcf_filtered = "{project_dir}data/vcf/filtered/{sample}_{lanes_normal}-{lanes_tumour}_somatic_filtered.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=".".join([get_lane(SAMPLES[SAMPLE]['G'][0]), get_lane(SAMPLES[SAMPLE]['G'][1])]), lanes_tumour=".".join([get_lane(SAMPLES[SAMPLE]['D'][0]), get_lane(SAMPLES[SAMPLE]['D'][1])]))
                if not vcf_somatic in VCF_SOMATIC:
                    VCF_SOMATIC.append(vcf_somatic)
                if not vcf_filtered in VCF_FILERED:
                    VCF_FILERED.append(vcf_filtered)
            elif (len(LANES)==3):
                file_3 = "{project_dir}data/bam/{sample}_{type}.{lane}.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane=get_lane(LANES[2]))
                if (os.path.isfile(file_3)):
                    MERGE_BAM.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bam".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    MERGE_BAI.append("{project_dir}data/bam/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge.bai".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    FASTQC.append("{project_dir}{fastq_dir}{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_fastqc.html".format(fastq_dir=config["FASTQC"]["DIR"], project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    VCF.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_for_pon.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    #VCF_IDX.append("{project_dir}data/vcf/{sample}_{type}.{lane_1}.{lane_2}.{lane_3}_marked_duplicates_BQSR_merge_for_pon.vcf.idx".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, type=TYPE, lane_1=get_lane(LANES[0]), lane_2=get_lane(LANES[1]), lane_3=get_lane(LANES[2])))
                    vcf_somatic = "{project_dir}data/vcf/{sample}_{lanes_normal}-{lanes_tumour}_somatic.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=".".join([get_lane(SAMPLES[SAMPLE]['G'][0]), get_lane(SAMPLES[SAMPLE]['G'][1]), get_lane(SAMPLES[SAMPLE]['G'][2])]), lanes_tumour=".".join([get_lane(SAMPLES[SAMPLE]['D'][0]), get_lane(SAMPLES[SAMPLE]['D'][1]), get_lane(SAMPLES[SAMPLE]['D'][2])]))
                    vcf_filtered = "{project_dir}data/vcf/filtered/{sample}_{lanes_normal}-{lanes_tumour}_somatic_filtered.vcf.gz".format(project_dir=config["PROJECT_DIR"], sample=SAMPLE, lanes_normal=".".join([get_lane(SAMPLES[SAMPLE]['G'][0]), get_lane(SAMPLES[SAMPLE]['G'][1]), get_lane(SAMPLES[SAMPLE]['G'][2])]), lanes_tumour=".".join([get_lane(SAMPLES[SAMPLE]['D'][0]), get_lane(SAMPLES[SAMPLE]['D'][1]), get_lane(SAMPLES[SAMPLE]['D'][2])]))
                    if not vcf_somatic in VCF_SOMATIC:
                        VCF_SOMATIC.append(vcf_somatic)
                    if not vcf_filtered in VCF_FILERED:
                        VCF_FILERED.append(vcf_filtered)

#print(MERGE)
TARGETS.extend(MERGE_BAM)
TARGETS.extend(MERGE_BAI)
TARGETS.extend(FASTQC)
TARGETS.extend(VCF_SOMATIC)
TARGETS.extend(VCF_FILERED)
#TARGETS.extend([config["PON_VCF"]])
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
        config["PROJECT_DIR"] + config["FASTQC"]["DIR"] + "{merged_samples}_marked_duplicates_BQSR_merge_fastqc.html"
    params:
        dir=config["FASTQC"]["DIR"],
        name="fastq_{merged_samples}",
        nthread=4
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o {params.dir}"


rule variant_calling_Mutect2:
    input:
        normal_bam = config["PROJECT_DIR"] + "data/bam/{sample}_G.{lanes_normal}_marked_duplicates_BQSR_merge.bam",
        tumour_bam = config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes_tumour}_marked_duplicates_BQSR_merge.bam",
	normal_bai = config["PROJECT_DIR"] + "data/bam/{sample}_G.{lanes_normal}_marked_duplicates_BQSR_merge.bai",
	tumour_bai = config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes_tumour}_marked_duplicates_BQSR_merge.bai"
    output:
        vcf_gz = config["PROJECT_DIR"] + "data/vcf/{sample}_{lanes_normal}-{lanes_tumour}_somatic.vcf.gz",
        f1r2_gz = config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_f1r2.tar.gz"
    wildcard_constraints:
        lanes_tumour="[0-9]\.[0-9]",
        lanes_normal="[0-9]\.[0-9]",
    params:
        ref=config["reference_GRCh37-lite"],
        PON=config["PON_VCF"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["intervals_list"],
        name="Mutect2_somatic_{sample}",
        nthread=config["mutect2"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk Mutect2 \
        -R {params.ref} \
        -L {params.intervals} \
        -I {input.normal_bam} \
        -I {input.tumour_bam} \
        -normal {wildcards.sample}_G_FREQEXCAP \
        --germline-resource {params.gnomad} \
        --panel-of-normals {params.PON} \
        --f1r2-tar-gz {output.f1r2_gz} \
        -O {output.vcf_gz}"


rule Mutect2_tumour_only:
    input:
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bam",
        bai=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bai"
    output:
        temp(config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.gz")
    params:
        ref=config["reference_GRCh37-lite"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
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
        -V gendb://{params.db} \
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


rule Calculate_Contamination_GetPileupSummaries:
    input:
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lanes}_marked_duplicates_BQSR_merge.bam",
	bai=config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lanes}_marked_duplicates_BQSR_merge.bai"
    output:
        config["PROJECT_DIR"] + "data/pileups/{sample}_{type}.{lanes}_pileups.table"
    wildcard_constraints:
        lanes="[0-9]\.[0-9]",
    params:
        gnomad=config["mutect2"]["gnomad"]["files"]["biallelic"],
	intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}_{type}.{lanes}",
        nthread=10
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.gnomad} \
            -L {params.intervals} \
            -O {output}"


rule Calculate_Contamination:
    input:
        tumour=config["PROJECT_DIR"] + "data/pileups/{sample}_D.{tumour_lanes}_pileups.table",
        matched=config["PROJECT_DIR"] + "data/pileups/{sample}_G.{normal_lanes}_pileups.table"
    output:
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}_{normal_lanes}-{tumour_lanes}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}_{normal_lanes}-{tumour_lanes}.tumour_segmentation.tsv"
    params:
        name="CalculateContamination_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk CalculateContamination \
            -I {input.tumour} \
            -matched {input.matched} \
            --tumor-segmentation {output.segmentation} \
            -O {output.contamination_table}"

rule LearnReadOrientationModel:
    input:
        config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_f1r2.tar.gz"
    output:
        config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_read-orientation-model.tar.gz"
    params:
        name="LearnReadOrientationModel_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk LearnReadOrientationModel \
            -I {input} \
            -O {output}"


rule GetPileupSummaries:
    input:
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes}_marked_duplicates_BQSR_merge.bam",
	bai=config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes}_marked_duplicates_BQSR_merge.bai"
    output:
        config["PROJECT_DIR"] + "data/f1r2/pileups/{sample}_D.{lanes}_getpileupsummaries.table"
    params:
        gnomad=config["mutect2"]["gnomad"]["files"]["biallelic"],
	intervals=config["intervals_list"], 
        name="GetPileupSummaries_{sample}",
        nthread=5
    conda: "../envs/gatk4.yaml"
    shell:
        "gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.gnomad}  \
            -L {params.intervals} \
            -O {output}"


rule FilterMutectCalls:
    input:
        vcf=config["PROJECT_DIR"] + "data/vcf/{sample}_{normal_lanes}-{tumour_lanes}_somatic.vcf.gz",
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}_{normal_lanes}-{tumour_lanes}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}_{normal_lanes}-{tumour_lanes}.tumour_segmentation.tsv",
        orientation=config["PROJECT_DIR"] + "data/f1r2/{sample}_{normal_lanes}-{tumour_lanes}_read-orientation-model.tar.gz"
    output:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered.vcf.gz"
    params:
        reference=config["reference_GRCh37-lite"],
        name="FilterMutectCalls_{sample}",
        nthread=config["FilterMutectCalls"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk FilterMutectCalls \
        -R {params.reference} \
        -V {input.vcf} \
        --contamination-table {input.contamination_table} \
        --tumor-segmentation {input.segmentation} \
        --orientation-bias-artifact-priors {input.orientation} \
        -O {output}"
