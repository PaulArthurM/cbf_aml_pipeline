
rule somatic_sniper:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam"
    output:
        "results/variantCalling/Somatic-sniper/{sample}_somatic-sniper.snv"
    params:
        name="somatic-sniper_{sample}",
        nthread=5,
        ref = config["reference_GRCh37-lite"]
    conda:
        "../envs/somaticSniper.yaml"
    shell:
        "bam-somaticsniper \
            -q 40 \
            -Q 40 \
            -f {params.ref} \
            {input.tumor} \
            {input.normal} \
            {output}"
