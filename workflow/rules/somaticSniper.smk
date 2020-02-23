
rule somatic_sniper:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam"
    output:
        "results/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.snv"
    params:
        name="somatic-sniper_{sample}",
        nthread=1,
        ref = config["reference_GRCh37-lite"]
    conda:
        "../envs/somaticSniper.yaml"
    shell:
        "bam-somaticsniper \
            -F classic \
            -f {params.ref} \
            {input.tumor} \
            {input.normal} \
            {output}"
