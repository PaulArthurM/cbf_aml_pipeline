
rule VariantFiltration:
    input:
        "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered_twice.vcf.gz"
    params:
        name="VariantFiltration_{sample}",
        nthread=5l,
        ref=config["reference_GRCh37-lite"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        'gatk VariantFiltration \
            --output {output} \
            --variant {input} \
            --filterName VariantAlleleCount    --filterExpression "VariantAlleleCount < 3" \
            --filterName VariantCountControl   --filterExpression "VariantAlleleCountControl > 1" \
            --filterName VariantBaseQualMedian --filterExpression "VariantBaseQualMedian < 25.0" \
            --filterName VariantMapQualMedian  --filterExpression "VariantMapQualMedian < 40.0" \
            --filterName MapQualDiffMedian     --filterExpression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
            --filterName LowMapQual            --filterExpression "LowMapQual > 0.05"'
