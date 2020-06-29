rule somatic_sniper:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam"
    output:
        "results/{token}/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf"
    params:
        name="somatic-sniper_{sample}",
        nthread=5,
        ref = config["reference"],
        #extra=extra_somatic_sniper
    log:
        "logs/{token}/somatic_sniper/{sample}.log"
    conda:
        "../envs/somaticSniper.yaml"
    shell:
        "bam-somaticsniper \
            -F vcf \
            -f {params.ref} \
            -q 20 \
            -Q 20 \
            {input.tumor} \
            {input.normal} \
            {output}"
