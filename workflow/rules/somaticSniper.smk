
rule somatic_sniper:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam"
    output:
        "results/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf"
    params:
        name="somatic-sniper_{sample}",
        nthread=1,
        ref = config["reference_GRCh37-lite"]
    conda:
        "../envs/somaticSniper.yaml"
    shell:
        "bam-somaticsniper \
            -F vcf \
            -q 20 \
            -Q 20 \
            -f {params.ref} \
            {input.tumor} \
            {input.normal} \
            {output}"
