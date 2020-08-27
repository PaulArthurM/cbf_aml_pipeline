# A Snakefile for a Variant Calling Pipeline integrating multiple variant callers.



import re
import numpy as np
import pandas as pd

# Load yaml configuration file
configfile: "config/config_lobry_serveur.yaml"


# Load PEP files
pepfile: "pep/config.yaml"


# Include all rules's scripts
include: "workflow/rules/commons.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/utils.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/mutect2.smk"
include: "workflow/rules/panelsOfNormals.smk"
include: "workflow/rules/strelka.smk"
include: "workflow/rules/freebayes.smk"
include: "workflow/rules/somaticSniper.smk"
include: "workflow/rules/filteringVCF.smk"
include: "workflow/rules/varlociraptor.smk"
include: "workflow/rules/varscan.smk"
include: "workflow/rules/sequenza.smk"



rule all:
    # See get_input function in workflow/rules/commons.smk
    # Define all output for Snakemake
    input: get_input
