# A Snakefile for a pre-processing of paired normal/tumour samples.


import json
import re
import os.path
from os import listdir
from os.path import isfile, join

include: "rules/functions.smk"
include: "rules/preprocessing.smk"
include: "rules/variantCalling.smk"
include: "rules/panelsOfNormals.smk"
include: "rules/annotation.smk"
include: "rules/utils.smk"

# Configuration file
configfile: "config/config.yaml"

# Load json configuration file
CONFIG_JSON = json.load(open(config["SAMPLES"]))

SAMPLES = CONFIG_JSON['samples']


rule all:
    # See get_input function in workflow/rules/functions.smk
    input: get_input
