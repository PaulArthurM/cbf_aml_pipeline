# A Snakefile for a Variant Calling Pipeline integrating multiple variant callers.

from snakemake.utils import min_version
min_version("5.24.0")


import re
import numpy as np
import pandas as pd

# Load yaml configuration file
configfile: "config/config_lobry_serveur.yaml"


# Load PEP files
pepfile: "pep/config.yaml"


# Include all rules's scripts
include: "workflow/rules/commons.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/utils.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/mutect2.smk"
include: "workflow/rules/panel_of_normals.smk"
include: "workflow/rules/strelka2.smk"
include: "workflow/rules/freebayes.smk"
include: "workflow/rules/somatic_sniper.smk"
include: "workflow/rules/filtering.smk"
include: "workflow/rules/varscan.smk"
include: "workflow/rules/sequenza.smk"



rule all:
    # See get_input function in workflow/rules/commons.smk
    # Define all output for Snakemake
    input:
        get_input
