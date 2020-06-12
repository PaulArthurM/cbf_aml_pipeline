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
    SAMPLES = sample_sheet['samples']
    LANES = SAMPLES[wildcards.sample][wildcards.type]
    fileToMerge = ""
    for file in getBamToMerge(wildcards):
        fileToMerge += " -I " + str(file)
    return fileToMerge


def getBamToMerge(wildcards):
    """Return a list containing all files expected as input for MergeSamFiles command."""
    SAMPLES = sample_sheet['samples']
    out = []
    SAMPLE = sample_sheet[sample_sheet.samples.eq(wildcards.sample)]
    if wildcards.type == "D":
        for bam in SAMPLE.at[0, "somatic_path"].split(" "):
            template = "results/preprocessing/" + wildcards.sample + "_" + wildcards.type + "." + getLane(bam) + "_marked_duplicates_BQSR.bam"
            out.append(template)
    if wildcards.type == "G":
        for bam in SAMPLE.at[0, "germline_path"].split(" "):
            template = "results/preprocessing/" + wildcards.sample + "_" + wildcards.type + "." + getLane(bam) + "_marked_duplicates_BQSR.bam"
            out.append(template)
    return out


def get_sample_name(sample):
    """Return sample corresponding to given file."""
    return re.match("(.+?)\.bam$", sample).group(1)


def get_id(sample):
    return re.search("-(\w+)\.", sample).group(1)


def get_lane(sample):
    return re.match("\.(\d+)$", sample).group(1)


def getLane(sample):
    return re.search("\.(\d)\.bam$", sample).group(1)


def get_tumour_lanes(sample):
    return ".".join([get_lane(lane) for lane in sample['D']])


def get_normal_lanes(sample):
    return ".".join([get_lane(lane) for lane in sample['G']])

def get_input(wildcards):
    wanted_input = []
    # Load json configuration file
    SAMPLES = sample_sheet['samples']
    TOKEN = config["token"]
    if config["panelsOfNormals"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/pon/{sample}_G_marked_duplicates_BQSR_merge_for_pon.vcf.gz", sample=SAMPLES, type=['G', 'D'], token=TOKEN))
    if config["mutect2"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz", sample=SAMPLES, token=TOKEN))
    if config["strelka"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/strelka/{sample}/strelka_calls.vcf.gz", sample=SAMPLES, token=TOKEN))
    if config["freebayes"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/freebayes/{sample}/freebayes_calls.vcf", sample=SAMPLES, token=TOKEN))
    if config["somaticSniper"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf", sample=SAMPLES, token=TOKEN))
    if config["annovar"]["activate"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/annovar/{sample}.hg19_multianno.vcf", sample=SAMPLES, token=TOKEN))
    if config["FASTQC"]["activate"] == True:
        #wanted_input.extend(expand("results/{token}/quality_control/{sample}_{type}_fastqc.html", sample=SAMPLES, type=['G', 'D'], token=TOKEN))
        wanted_input.append("results/{token}/report/multiqc_report.html".format(token=TOKEN))
    #if config['sequenza']['activate'] == True:
        #wanted_input.extend(expand('results/{token}/sequenza/{sample}.small.seqz.gz', sample=SAMPLES, token=TOKEN))
        #wanted_input.extend(expand('results/{token}/sequenza/{sample}_seqz/{sample}_segments.bed', sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["diploid_variants"] == True:
        wanted_input.extend(expand('results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass_diploid.vcf', sample=SAMPLES, token=TOKEN))
    if config["VariantFiltering"]["oxog_filtering"] == True:
        wanted_input.extend(expand("results/{token}/variantCalling/vcf/mutect2/oxog_filtered/{sample}_oxog_filtered.vcf.gz", sample=SAMPLES, token=TOKEN))
    return wanted_input
