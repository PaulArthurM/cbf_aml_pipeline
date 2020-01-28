rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bam",
        normal_index = config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bai",
        tumor = config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bam",
        tumor_index = config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bai",
    output:
        directory(config["PROJECT_DIR"] + "results/variantCalling/strelka/raw/{sample}_strelka_vcf")
    params:
        name="strelka_{sample}",
        nthread=5,
        ref = config["reference_GRCh37-lite"]
    conda:
        "workflow/envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir {output} \
            && \
            {output}/runWorkflow.py \
            --mode local \
            --jobs {params.nthread}"
