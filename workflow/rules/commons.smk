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
    for bam in SAMPLES[wildcards.sample][wildcards.type]:
        template = "results/preprocessing/" + wildcards.sample + "_" + wildcards.type + "." + get_lane(bam) + "_marked_duplicates_BQSR.bam"
        out.append(template)
    return out


def get_sample_name(sample):
    """Return sample corresponding to given file."""
    return re.match("(.+?)\.bam$", sample).group(1)


def get_id(sample):
    return re.search("-(\w+)\.", sample).group(1)


def get_lane(sample):
    return re.search("\.(\d+)$", sample).group(1)


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
    if config["panelsOfNormals"]["activate"] == True:
        wanted_input.extend(expand("results/pon/{sample}_{type}_marked_duplicates_BQSR_merge_for_pon.vcf.gz", sample=SAMPLES, type=['G', 'D']))
    if config["mutect2"]["activate"] == True:
        wanted_input.extend(expand("results/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered_fdr05.vcf.gz", sample=SAMPLES))
        wanted_input.extend(expand("results/variantCalling/vcf/mutect2/oxog_filtered/{sample}_oxog_filtered.vcf.gz", sample=SAMPLES))
    if config["strelka"]["activate"] == True:
        wanted_input.extend(expand("results/variantCalling/strelka/{sample}/strelka_calls.vcf.gz", sample=SAMPLES))
    if config["freebayes"]["activate"] == True:
        wanted_input.extend(expand("results/variantCalling/freebayes/{sample}/freebayes_calls.vcf", sample=SAMPLES))
    if config["somaticSniper"]["activate"] == True:
        wanted_input.extend(expand("results/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf", sample=SAMPLES))
    if config["annovar"]["activate"] == True:
        wanted_input.extend(expand("results/variantCalling/annovar/{sample}.hg19_multianno.vcf", sample=SAMPLES))
    if config["FASTQC"]["activate"] == True:
        wanted_input.extend(expand("results/quality_control/{sample}_{type}_fastqc.html", sample=SAMPLES, type=['G', 'D']))
        wanted_input.append("results/report/multiqc_report.html")
    if config['sequenza']['activate'] == True:
        wanted_input.extend(expand('results/sequenza/small.{sample}.seqz.gz', sample=SAMPLES))
        wanted_input.extend(expand('results/sequenza/{sample}_seqz/{sample}_segments.bed', sample=SAMPLES))
    return wanted_input
