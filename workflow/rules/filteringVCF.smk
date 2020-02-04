
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
            --filter-name VariantAlleleCount    --filterExpression "VariantAlleleCount < 3" \
            --filter-name VariantCountControl   --filterExpression "VariantAlleleCountControl > 1" \
            --filter-name VariantBaseQualMedian --filterExpression "VariantBaseQualMedian < 25.0" \
            --filter-name VariantMapQualMedian  --filterExpression "VariantMapQualMedian < 40.0" \
            --filter-name MapQualDiffMedian     --filterExpression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
            --filter-name LowMapQual            --filterExpression "LowMapQual > 0.05"'
