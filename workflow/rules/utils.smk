rule MergeSamFiles:
    input:
        getBamToMerge
    output:
        "results/preprocessing/{sample, [A-Za-z0-9]+}_{type}.bam"
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


rule IndexFeatureFile:
    input:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz.tbi"
    params:
        name="IndexFeatureFile_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk IndexFeatureFile -F {input}"
