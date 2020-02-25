

rule mpileup:
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/preprocessing/{sample}_{type}.mpileup"
    params:
        ref = config["reference_GRCh37-lite"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools mpileup -f {params.ref} {input} > {output}"


rule varscan:
    input:
        normal="results/preprocessing/{sample}_G.mpileup",
        tumor="results/preprocessing/{sample}_D.mpileup"
    output:
        "results/variantCalling/varscan/{sample}/varscan_calls.vcf"
    params:
        ""
    conda:
        "../envs/varscan.yaml"
    shell:
        "varscan somatic \
            {input.normal} \
            {input.tumor} \
            {output} \
            --output-vcf \
            --min-coverage 1 && \
            mv {output}.snp {output}"

"""
rule index_vcf:
    input:
        "results/variantCalling/{tool}/{sample}/{tool}_calls.vcf"
    output:
        "results/variantCalling/{tool}/{sample}/{tool}_calls.vcf.gz"  # either .vcf or .bcf
    shell:
        "bgzip {input}"
"""
