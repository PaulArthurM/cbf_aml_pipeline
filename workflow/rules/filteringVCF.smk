
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
            --filter-name VariantAlleleCount    --filter-expression "VariantAlleleCount < 3" \
            --filter-name VariantCountControl   --filter-expression "VariantAlleleCountControl > 1" \
            --filter-name VariantBaseQualMedian --filter-expression "VariantBaseQualMedian < 25.0" \
            --filter-name VariantMapQualMedian  --filter-expression "VariantMapQualMedian < 40.0" \
            --filter-name MapQualDiffMedian     --filter-expression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
            --filter-name LowMapQual            --filter-expression "LowMapQual > 0.05"'
