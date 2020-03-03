def extra_somatic_sniper(wildcards):
    configfile: "config/config.yaml"
    extra = config['somaticSniper']['extra']
    return extra


rule somatic_sniper:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam"
    output:
        "results/variantCalling/somatic-sniper/{sample}/somatic-sniper_calls.vcf"
    params:
        name="somatic-sniper_{sample}",
        nthread=5,
        ref = config["reference_GRCh37-lite"],
        extra=extra_somatic_sniper
    conda:
        "../envs/somaticSniper.yaml"
    shell:
        "bam-somaticsniper \
            -F vcf \
            {params.extra} \
            -f {params.ref} \
            {input.tumor} \
            {input.normal} \
            {output}"
