

"results/preprocessing/{sample}_G.bam"
"results/preprocessing/{sample}_D.bam"


rule preprocessing:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
        vcf="results/variantCalling/{tool}/{sample}/{tool}_calls.vcf.gz"
    output:
        "results/varlociraptor/{sample}/varlociraptor_{tool}_{type}_preprocessing.bcf"
    params:
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
        ""
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants tumor-normal \
            --purity 0.75 \
            --tumor {input.tumor} \
            --normal {input.normal} > {output}"
