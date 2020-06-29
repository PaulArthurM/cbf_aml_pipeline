rule varlociraptor_preprocessing:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
        vcf="results/{token}/variantCalling/{tool}/{sample}/{tool}_calls.vcf.gz"
    output:
        "results/{token}/varlociraptor/{sample}/varlociraptor_{tool}_{type}_preprocessing.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}_{tool}",
        nthread=5,
        ref=config["reference"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.ref} \
        --bam {input.bam} \
        --output {output} < {input.vcf}"



rule varlociraptor_calling:
    input:
        tumor="results/{token}/varlociraptor/{sample}/varlociraptor_{tool}_D_preprocessing.bcf",
        normal="results/{token}/varlociraptor/{sample}/varlociraptor_{tool}_G_preprocessing.bcf"
    output:
        "results/{token}/varlociraptor/{sample}/{tool}_calls.bcf"
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



rule varlociraptor_filter:
    input:
        "results/{token}/varlociraptor/{sample}/{tool}_calls.bcf"
    output:
        "results/{token}/varlociraptor/{sample}/{tool}_calls.filtered.bcf"
    params:
        name="varlociraptor_filter_{sample}_{tool}",
        nthread=5
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr \
        --events SOMATIC_TUMOR \
        --fdr 0.05 \
        --var SNV {input} > {output}"
