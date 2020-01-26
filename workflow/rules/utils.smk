
def getBamToMergeCommand(wildcards):
    # Configuration file
    configfile: "config/config.yaml"
    CONFIG_JSON = json.load(open(config["SAMPLES"]))
    SAMPLES = CONFIG_JSON['samples']
    LANES = SAMPLES[wildcards.sample][wildcards.type]
    lanesToMerge = ""
    for lane in LANES:
        lanesToMerge += "-I " + str(lane)
    return lanesToMerge


def getBamToMerge(wildcards):
    configfile: "config/config.yaml"
    CONFIG_JSON = json.load(open(config["SAMPLES"]))
    SAMPLES = CONFIG_JSON['samples']
    return SAMPLES[wildcards.sample][wildcards.type]


def test(wildcards):
    SAMPLES = CONFIG_JSON['samples']
    out = []
    template = config["PROJECT_DIR"] + "data/bam/" + wildcards.sample + "_" + wildcards.type + "_{lane}_marked_duplicates_BQSR.bam"#.format(sample=wildcards.sample, type=wildcards.type)
    lanes = [getLane(bam) for bam in SAMPLES[wildcards.sample][wildcards.type]]
    out.extend(expand(template, lane=lanes))
    return out



rule merge_bam:
    input:
        test#getBamToMerge
    output:
        config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bam"
    conda:
        "../envs/gatk4.yaml"
    params:
        name="merge_{sample}",
        nthread=5,
        bamToMerge=getBamToMergeCommand
    shell:
        "gatk MergeSamFiles \
            {params.bamToMerge} \
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
        config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bam"
    output:
        config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bai"
    conda:
        "../envs/samtools.yaml"
    params:
        name="index_{sample}_{type}",
        nthread=config["samtools"]["nthread"]
    shell:
        "samtools index -b {input} {output}"


rule fastqc:
    input:
        config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bam"
    output:
        config["PROJECT_DIR"] + config["FASTQC"]["DIR"] + "{sample}_{type}_fastqc.html"
    params:
        dir=config["FASTQC"]["DIR"],
        name="fastq_{sample}_{type}",
        nthread=4
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o {params.dir}"


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
