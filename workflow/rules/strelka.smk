
"""
# configuration
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam normal.bam \
    --tumorBam tumor.bam \
    --referenceFasta hg38.fa \
    --runDir demo_somatic
# execution on a single local machine with 20 parallel jobs
demo_somatic/runWorkflow.py -m local -j 20
"""

def get_strelka_input_normal_bam(wildcards):
    return config["PROJECT_DIR"] + "data/bam/{sample}_G.{lanes_normal}_marked_duplicates_BQSR_merge.bam".format(sample = wildcards.sample, lanes_normal = wildcards.lanes_normal)


def get_strelka_input_tumour_bam(wildcards):
    return config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes_tumour}_marked_duplicates_BQSR_merge.bam".format(sample = wildcards.sample, lanes_tumour = wildcards.lanes_tumour)


rule configureStrelkaSomaticWorkflow:
    input:
        normalBam = get_strelka_input_normal_bam,
        tumorBam = get_strelka_input_tumour_bam
    params:
        ref=config["reference_GRCh37-lite"],
        name="Mutect2_somatic_{sample}",
        pathToConfigure=["strelka"]["STRELKA_INSTALL_PATH"],
        runDir=["strelka"]["runDir"],
        nthread=config["mutect2"]["nthread"]
    output:

    shell:
        "${params.STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
            --normalBam {input.normalBam} \
            --tumorBam {input.tumorBam} \
            --referenceFasta {params.ref} \
            --runDir {params.runDir}"
