# A Snakefile for a pre-processing of paired normal/tumour samples.


import json
import re
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

# Load json configuration file
configfile: "config/config_laptop.yaml"


# Load tsv sample sheet
sample_sheet = pd.read_csv(config['sample_sheet'], sep=";")


include: "workflow/rules/commons.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/utils.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/mutect2.smk"
include: "workflow/rules/panelsOfNormals.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/strelka.smk"
include: "workflow/rules/freebayes.smk"
include: "workflow/rules/somaticSniper.smk"
include: "workflow/rules/filteringVCF.smk"
include: "workflow/rules/varlociraptor.smk"
include: "workflow/rules/varscan.smk"
include: "workflow/rules/sequenza.smk"


rule all:
    # See get_input function in workflow/rules/commons.smk
    input: get_input
