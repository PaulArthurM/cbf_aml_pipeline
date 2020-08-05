# This script contain rules for:
#    - merging bam files from one sample in one nex BAM file with MergeSamFiles
#    - indexing bam files in bai files with samtools index
#    - indexing vcf.gz files in vcf.gz.tbi files with IndexFeatureFile



rule MergeSamFiles:
    """Run MergeSamFiles on a list of bam files and return a merged bam file"""
    input:
        getBamToMerge
    output:
        "results/preprocessing/{sample, [A-Za-z0-9]+}_{type}.bam"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    params:
        name="merge_{sample}_{type}",
        nthread=5,
        bamToMerge=getBamToMergeCommand
    log:
        "logs/preprocessing/MergeSamFiles/{sample}_{type}.log"
    wildcard_constraints:
        type="[DG]"
    shell:
        "gatk MergeSamFiles \
            {params.bamToMerge} \
            -O {output} 2> {log}"



rule samtools_index:
    """Run samtools index on a bam file and return an index bai file"""
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/preprocessing/{sample}_{type}.bai"
    conda:
        "../envs/samtools.yaml"
    params:
        name="index_{sample}_{type}",
        nthread=config["samtools"]["nthread"]
    log:
        "logs/preprocessing/samtools_index/{sample}_{type}.log"
    shell:
        "samtools index {input} {output} 2> {log}"


rule IndexFeatureFile:
    """Run IndexFeatureFile on a vcf.gz file and return an index vcf.gz.tbi file"""
    input:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz.tbi"
    params:
        name="IndexFeatureFile_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk IndexFeatureFile -F {input}"



rule biallelic_vcf:
    input:
        config["mutect2"]["gnomad"]["files"]["raw"]
    output:
        "ressources/somatic-b37_af-only-gnomad.b37.BIALLELIC.vcf"
    params:
        name="SelectVariants_BIALLELIC",
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk SelectVariants \
        -R {params.ref} \
        -V {input} \
        -O {output} \
        --restrict-alleles-to BIALLELIC 2> {log}"
