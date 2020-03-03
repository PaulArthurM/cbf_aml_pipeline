# A Snakefile for a pre-processing of paired normal/tumour samples.


import json
import re
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

# Load json configuration file
configfile: "config/config.yaml"
CONFIG_JSON = json.load(open(config["SAMPLES"]))
SAMPLES = CONFIG_JSON['samples']

# Load tsv sample sheet
sample_sheet = pd.read_csv(config['sample_sheet'])

include: "workflow/rules/utils.smk"
include: "workflow/rules/functions.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/mutect2.smk"
#include: "workflow/rules/panelsOfNormals.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/strelka.smk"
include: "workflow/rules/freebayes.smk"
include: "workflow/rules/somaticSniper.smk"
include: "workflow/rules/filteringVCF.smk"
include: "workflow/rules/varlociraptor.smk"
include: "workflow/rules/varscan.smk"


def get_input(wildcards):
    wanted_input = []
    # Load json configuration file
    CONFIG_JSON = json.load(open(config["SAMPLES"]))
    #SAMPLES = CONFIG_JSON['samples']
    SAMPLES = sample_sheet['samples']
    wanted_input.extend(expand("results/preprocessing/{sample}_{type}.bam", sample=SAMPLES, type=['G', 'D']))
    wanted_input.extend(expand("results/preprocessing/{sample}_{type}.bai", sample=SAMPLES, type=['G', 'D']))
    if config["panelsOfNormals"]["to_use"] == True:
        wanted_input.extend(expand("results/pon/{sample}_{type}_marked_duplicates_BQSR_merge_for_pon.vcf.gz", sample=SAMPLES, type=['G', 'D']))
    if config["mutect2"]["to_use"] == True:
        #wanted_input.extend(expand("results/variantCalling/mutect2/{sample}/mutect2_calls.vcf.gz", sample=SAMPLES))
        #wanted_input.extend(expand("results/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz", sample=SAMPLES))
        wanted_input.extend(expand("results/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf", sample=SAMPLES))
    if config["strelka"]["to_use"] == True:
        wanted_input.extend(expand("results/variantCalling/Strelka/{sample}/results/variants/somatic.snvs.vcf.gz", sample=SAMPLES))
    if config["freebayes"]["to_use"] == True:
        wanted_input.extend(expand("results/variantCalling/freebayes/raw/{sample}_freebayes.vcf", sample=SAMPLES))
    if config["somaticSniper"]["to_use"] == True:
        wanted_input.extend(expand("results/variantCalling/Somatic-sniper/{sample}_somatic-sniper.snv", sample=SAMPLES))
    if config["annovar"]["to_use"] == True:
        wanted_input.extend(expand("results/variantCalling/annovar/{sample}.hg19_multianno.vcf", sample=SAMPLES))
    if config["FASTQC"]["to_use"] == True:
        wanted_input.extend(expand(config["FASTQC"]["DIR"] + "{sample}_{type}_fastqc.html", sample=SAMPLES, type=['G', 'D']))
    return wanted_input


rule all:
    # See get_input function in workflow/rules/functions.smk
    input: get_input
