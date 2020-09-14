# This script contain functions for:
#    - extra parameters for various tools
#    - returns command's parameters for MergeSamFiles
#    - returns list of input files for MergeSamFiles
#    - extract some informations from sample names with regex
#    - and also the main function allowing to define input of rule 'all' thanks
#      to config file and sample sheet.



def extra_freebayes(wildcards):
    """Returns a string containing additional parameters for the Freebayes command."""
    configfile: "config/config.yaml"
    extra = config['freebayes']['extra']
    return extra


def extra_mutect2(wildcards):
    """Returns a string containing additional parameters for the Mutect2 command."""
    configfile: "config/config.yaml"
    extra = config['mutect2']['extra']
    return extra


def extra_somatic_sniper(wildcards):
    """Returns a string containing additional parameters for the Somatic Sniper command."""
    configfile: "config/config.yaml"
    extra = config['somaticSniper']['extra']
    return extra


def extra_strelka(wildcards):
    """Returns a string containing additional parameters for the Strelka command."""
    configfile: "config/config.yaml"
    extra = config['strelka']['extra']
    return extra


def getBamToMergeCommand(wildcards):
    """Return a string to inject in MergeSamFiles command for multiple input files."""
    SAMPLES = pep.sample_table["sample_name"]
    fileToMerge = ""
    for file in getBamToMerge(wildcards):
        fileToMerge += " -I " + str(file)
    return fileToMerge


# def getBamToMerge(wildcards):
#     """Return a list containing all files expected as input for MergeSamFiles command."""
#     SAMPLES = pep.sample_table["sample_name"]
#     out = []
#     SAMPLE = sample_sheet.set_index("samples", drop = False)
#     if wildcards.type == "D":
#         for bam in SAMPLE.at[wildcards.sample, "somatic_path"].split(" "):
#             template = "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates_BQSR.bam".format(sample=wildcards.sample, type=wildcards.type, lane=getLane(bam))
#             out.append(template)
#     if wildcards.type == "G":
#         for bam in SAMPLE.at[wildcards.sample, "germline_path"].split(" "):
#             template = "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates_BQSR.bam".format(sample=wildcards.sample, type=wildcards.type, lane=getLane(bam))
#             out.append(template)
#     return out


def getBamToMerge(wildcards):
    """Return a list containing all files expected as input for MergeSamFiles command."""
    SAMPLE_TABLE = pep.sample_table["sample_name"]
    SUBSAMPLES_TABLE = pep.subsample_table
    out = []
    if wildcards.type == "D":
        for lane in SUBSAMPLES_TABLE[(SUBSAMPLES_TABLE.sample_name == wildcards.sample) & (SUBSAMPLES_TABLE.type != "D")]["lane"]:
            template = "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates_BQSR.bam".format(sample=wildcards.sample, type=wildcards.type, lane=lane)
            out.append(template)
    if wildcards.type == "G":
        for lane in SUBSAMPLES_TABLE[(SUBSAMPLES_TABLE.sample_name == wildcards.sample) & (SUBSAMPLES_TABLE.type != "G")]["lane"]:
            template = "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates_BQSR.bam".format(sample=wildcards.sample, type=wildcards.type, lane=lane)
            out.append(template)
    return out


def get_sample_name(sample):
    """Return sample corresponding to given file."""
    return re.match("(.+?)\.bam$", sample).group(1)


def get_id(sample):
    """Return sample's ID corresponding to given file."""
    return re.search("-(\w+)\.", sample).group(1)


# useless?
def get_lane(sample):
    """Return sample's lane corresponding to given file."""
    return re.match("\.(\d+)$", sample).group(1)


def getLane(sample):
    """Return sample's lane corresponding to given file."""
    return re.search("\.(\d)\.bam$", sample).group(1)


def get_tumour_lanes(sample):
    """Return tumor's lanes corresponding to given file."""
    return ".".join([get_lane(lane) for lane in sample['D']])


def get_normal_lanes(sample):
    """Return normal's lanes corresponding to given file."""
    return ".".join([get_lane(lane) for lane in sample['G']])


def get_input(wildcards):
    """Return a list of all input file based on config file and sample sheet."""
    wanted_input = []
    # Load sample sheet
    SAMPLES = pep.sample_table["sample_name"]
    #  Token is used to run the pipeline multiple times without erasing other results,
    #  just by changing token's value in config file.
    TOKEN = config["token"]

    if config["mode"] == "preprocessing":
        return expand("results/preprocessing/{sample}_{type}.bam", sample=SAMPLES, type=['G', 'D'])
    if config["mode"] == "qc":
        return "results/{token}/quality_control/report/multiqc_report.html".format(token=TOKEN)
    if config["mode"] == "variantCalling":
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz", sample=SAMPLES, token=TOKEN))
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/strelka/merged/{sample}.merged.vcf.gz", sample=SAMPLES, token=TOKEN))
        #wanted_input.extend(expand("results/{token}/variantCalling/freebayes/{sample}/freebayes_calls.vcf", sample=SAMPLES, token=TOKEN))
        return wanted_input


    if config["panelsOfNormals"]["activate"] == True:  #  should not be an option?
        wanted_input.extend(expand("results/{token}/pon/{sample}_G_marked_duplicates_BQSR_merge_for_pon.vcf.gz", sample=SAMPLES, type=['G', 'D'], token=TOKEN))
    if config["mutect2"]["activate"] == True:  # useless if vcf output from mutect2 is activated in VariantFiltering
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz", sample=SAMPLES, token=TOKEN))
    if config["strelka"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/strelka/merged/{sample}.merged.vcf.gz", sample=SAMPLES, token=TOKEN))
    if config["freebayes"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/freebayes/{sample}/freebayes_calls.vcf", sample=SAMPLES, token=TOKEN))
    if config["somaticSniper"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf", sample=SAMPLES, token=TOKEN))
    if config["FASTQC"]["activate"] == True:
        wanted_input.append("results/{token}/quality_control/report/multiqc_report.html".format(token=TOKEN))
    if config["VariantFiltering"]["pass"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf", sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["diploid_variants"] == True:
        wanted_input.extend(expand('results/{token}/variantCalling/vcf/mutect2/diploid/{sample}_somatic_filtered_diploid.vcf', sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["mutect2_strelka2_intersection"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/intersection/{sample}_mutect2_strelka2_intersection.vcf", sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["varlociraptor"]["classic"] == True:
        wanted_input.extend(expand("results/{token}/varlociraptor/{sample}/classic_calls.filtered.bcf", sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["varlociraptor"]["scenario"] == True:
        wanted_input.extend(expand("results/{token}/varlociraptor/{sample}/scenario_calls.filtered.bcf", sample=SAMPLES, token=TOKEN))
    if config["sequenza"]["activate"] == True:
        wanted_input.extend(expand('results/{token}/sequenza/{sample}_seqz/{sample}_segments.bed', sample=SAMPLES, token=TOKEN))

    return wanted_input
