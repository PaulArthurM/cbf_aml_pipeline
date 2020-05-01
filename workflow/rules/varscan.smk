

rule samtools_mpileup:
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/preprocessing/{sample}_{type}.mpileup"
    params:
        ref = config["reference"]
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools mpileup -f {params.ref} {input} > {output}"


rule varscan:
    input:
        normal="results/preprocessing/{sample}_G.mpileup",
        tumor="results/preprocessing/{sample}_D.mpileup"
    output:
        "results/{token}/variantCalling/varscan/{sample}/varscan_calls.vcf"
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
