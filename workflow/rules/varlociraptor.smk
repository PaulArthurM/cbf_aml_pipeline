rule varlociraptor_preprocessing:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
        vcf="results/{token}/variantCalling/intersection/{sample}_mutect2_strelka2_intersection.vcf"
        #vcf="results/{token}/variantCalling/{tool}/{sample}/{tool}_calls.vcf.gz"
    output:
        "results/{token}/varlociraptor/{sample}/varlociraptor_{type}_preprocessing.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}",
        nthread=5,
        ref=config["reference"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.ref} \
        --bam {input.bam} \
        --candidates {input.vcf} \
        --output {output}"



rule varlociraptor_calling_scenario:
    input:
        tumor="results/{token}/varlociraptor/{sample}/varlociraptor_D_preprocessing.bcf",
        normal="results/{token}/varlociraptor/{sample}/varlociraptor_G_preprocessing.bcf"
    output:
        "results/{token}/varlociraptor/{sample}/scenario_calls.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}",
        nthread=5,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants generic \
            --scenario config/varlociraptor_events.yaml \
            --obs tumor={input.tumor} normal={input.normal} \
            > {output}"



rule varlociraptor_calling_classic:
    input:
        tumor="results/{token}/varlociraptor/{sample}/varlociraptor_D_preprocessing.bcf",
        normal="results/{token}/varlociraptor/{sample}/varlociraptor_G_preprocessing.bcf"
    output:
        "results/{token}/varlociraptor/{sample}/classic_calls.bcf"
    params:
        name="varlociraptor_preprocessing_{sample}",
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
        "results/{token}/varlociraptor/{sample}/{type}_calls.bcf"
    output:
        "results/{token}/varlociraptor/{sample}/{type}_calls.filtered.bcf"
    params:
        name="varlociraptor_filter_{sample}",
        nthread=5
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr \
        --events SOMATIC_TUMOR \
        --fdr 0.05 \
        --var SNV {input} > {output}"
