rule preprocessing:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
        vcf="results/variantCalling/{tool}/{sample}/{tool}_calls.vcf.gz"
    output:
        "results/varlociraptor/{sample}/varlociraptor_{tool}_{type}_preprocessing.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}_{tool}",
        nthread=5,
        ref=config["reference_GRCh37-lite"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.ref} \
        --bam {input.bam} \
        --output {output} < {input.vcf}"


rule calling:
    input:
        tumor="results/varlociraptor/{sample}/varlociraptor_{tool}_D_preprocessing.bcf",
        normal="results/varlociraptor/{sample}/varlociraptor_{tool}_G_preprocessing.bcf"
    output:
        "results/varlociraptor/{sample}/{tool}_calls.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}_{tool}",
        nthread=5,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants tumor-normal \
            --purity 0.75 \
            --tumor {input.tumor} \
            --normal {input.normal} > {output}"

"""
rule ffpe:
    input:
    output:
    params:
    conda:
    shell:
        "varlociraptor call variants generic --scenario scenario.yaml --obs relapse=relapse.bcf tumor=tumor.bcf normal=normal.bam > calls.bcf"

rule filter:
    input:
    output:
    params:
    conda:
    shell:
        "varlociraptor filter-calls control-fdr calls.bcf --events SOMATIC_TUMOR --fdr 0.05 --var SNV"
"""
