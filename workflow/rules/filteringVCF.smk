
rule VariantFiltration:
    input:
        "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered_twice.vcf.gz"
    params:
        name="VariantFiltration_{sample}",
        nthread=5,
        ref=config["reference_GRCh37-lite"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        'gatk VariantFiltration \
            --output {output} \
            --variant {input} \
            --filter-name VariantOrientation    --filter-expression "ALT_F1R2 > 0" \
            --filter-name VariantOrientation    --filter-expression "ALT_F2R1 > 0" \
            --filter-name VariantBaseQualMedian --filter-expression "MBQ < 25.0" \
            --filter-name VariantMapQualMedian  --filter-expression "MMQ < 40.0" \
            --filter-name VariantTLOD           --filter-expression "TLOD > 8.5"'
