def getBamToMergeCommand(wildcards):
    SAMPLES = CONFIG_JSON['samples']
    LANES = SAMPLES[wildcards.sample][wildcards.type]
    fileToMerge = ""
    for file in getBamToMerge(wildcards):
        fileToMerge += " -I " + str(file)
    return fileToMerge



def getBamToMerge(wildcards):
    SAMPLES = CONFIG_JSON['samples']
    out = []
    for bam in SAMPLES[wildcards.sample][wildcards.type]:
        template = "results/preprocessing/" + wildcards.sample + "_" + wildcards.type + "." + get_lane(bam) + "_marked_duplicates_BQSR.bam"#.format(sample=wildcards.sample, type=wildcards.type)
        out.append(template)
    return out


rule merge_bam:
    input:
        getBamToMerge
    output:
        "results/preprocessing/{sample}_{type}.bam"
    conda:
        "../envs/gatk4.yaml"
    params:
        name="merge_{sample}_{type}",
        nthread=5,
        bamToMerge=getBamToMergeCommand
    wildcard_constraints:
        type="[DG]"
    shell:
        "gatk MergeSamFiles \
            {params.bamToMerge} \
            -O {output}"



# Rule for create index from BAM file with samtools index
rule samtools_index:
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/preprocessing/{sample}_{type}.bai"
    conda:
        "../envs/samtools.yaml"
    params:
        name="index_{sample}_{type}",
        nthread=config["samtools"]["nthread"]
    shell:
        "samtools index -b {input} {output}"


rule fastqc:
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        config["FASTQC"]["DIR"] + "{sample}_{type}_fastqc.html"
    params:
        dir=config["FASTQC"]["DIR"],
        name="fastq_{sample}_{type}",
        nthread=4
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o {params.dir}"

#
# rule unzip_gz:
#     input:
#         vcf_gz = config["PROJECT_DIR"] + "results/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.gz"
#     output:
#         vcf = config["PROJECT_DIR"] + "results/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf"
#     params:
#         name="gunzip_{sample}_G.{lane}",
#         nthread=1
#     shell:
#         "gunzip {input.vcf_gz}"
#
#
 rule IndexFeatureFile:
     input:
         "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
     output:
        "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz.idx"
     params:
         name="IndexFeatureFile_{sample}",
         nthread=5
     conda:
         "../envs/gatk4.yaml"
     shell:
         "gatk IndexFeatureFile -F {input}"
