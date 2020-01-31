rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = "results/preprocessing/{sample}_G.bam",
        normal_index = "results/preprocessing/{sample}_G.bai",
        tumor = "results/preprocessing/{sample}_D.bam",
        tumor_index = "results/preprocessing/{sample}_D.bai",
    output:
        directory("results/variantCalling/strelka/raw/{sample}_strelka_vcf")
    params:
        name="strelka_{sample}",
        nthread=5,
        ref = config["reference_GRCh37-lite"]
    conda:
        "../envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir {output} \
            && \
            {output}/runWorkflow.py \
            --mode local \
            --jobs {params.nthread}"
