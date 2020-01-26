def get_strelka_input_normal_bam(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bam".format(sample = wildcards.sample)


def get_strelka_input_tumour_bam(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bam".format(sample = wildcards.sample)


def get_strelka_input_normal_bai(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bai".format(sample = wildcards.sample)


def get_strelka_input_tumour_bai(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bai".format(sample = wildcards.sample)


rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = get_strelka_input_normal_bam,
        normal_index = get_strelka_input_normal_bai,
        tumor = get_strelka_input_tumour_bam,
        tumor_index = get_strelka_input_tumour_bai,
    output:
        directory(config["PROJECT_DIR"] + "results/variantCalling/strelka/raw/{sample}_strelka_vcf")
    params:
        name="strelka_{sample}",
        nthread=5,
        ref = config["reference_GRCh37-lite"]
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
