rule annotate:
    input:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered.vcf.gz"
    output:
        config["PROJECT_DIR"] + "data/vcf/annotated/{sample}_{normal_lanes}-{tumour_lanes}_variants.funcotated.vcf"
    params:
        reference=config["reference_GRCh37-lite"],
        data_sources=config["funcotator"]["directory"],
        name="Funcotator_annotate_{sample}",
        nthread=4
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk Funcotator \
        --variant {input} \
        --reference {params.reference} \
        --ref-version b37 \
        --data-sources-path {params.data_sources} \
        --output {output} \
        --allow-hg19-gencode-b37-contig-matching true \
        --output-file-format VCF"
