
rule FilterMutectCalls:
    input:
        vcf="results/{token}/variantCalling/mutect2/{sample}/mutect2_calls.vcf.gz",
        contamination_table="results/{token}/variantCalling/mutect2/pileups/contamination/{sample}.contamination.table",
        segmentation="results/{token}/variantCalling/mutect2/pileups/segmentation/{sample}.tumour_segmentation.tsv",
        orientation="results/{token}/variantCalling/mutect2/f1r2/{sample}_read-orientation-model.tar.gz"
    output:
        "results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    params:
        reference=config["reference"],
        name="FilterMutectCalls_{sample}",
        nthread=config["FilterMutectCalls"]["nthread"]
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk FilterMutectCalls \
        -R {params.reference} \
        -V {input.vcf} \
        --contamination-table {input.contamination_table} \
        --tumor-segmentation {input.segmentation} \
        --orientation-bias-artifact-priors {input.orientation} \
        --max-events-in-region 4 \
        -O {output}"



rule keep_pass_variants:
    input:
        "results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf"
    params:
        name="keep_pass_variants_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view \
        -f .,PASS \_GRCh37-lite \
        {input} > {output}"



rule bedtools_keep_diploid_variants:
    input:
        vcf="results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf",
        diploid_regions="results/{token}/sequenza/{sample}_seqz/{sample}_segments.bed"
    output:
        diploid_snv="results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass_diploid.vcf"
    params:
        name="Keep_diploid_{sample}",
        nthread=5,
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.vcf} -b {input.diploid_regions}"



rule VariantFiltration:
    input:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/{token}/variantCalling/mutect2/filtered/{sample}_somatic_filtered_twice.vcf.gz"
    params:
        name="VariantFiltration_{sample}",
        nthread=5,
        ref=config["reference"]
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        'gatk VariantFiltration \
            --output {output} \
            --variant {input} \
            --filter-name VariantOrientation    --filter-expression "F1R2 < 4.0" \
            --filter-name VariantOrientation    --filter-expression "F2R1 < 4.0" \
            --filter-name VariantBaseQualMedian --filter-expression "MBQ < 25.0" \
            --filter-name VariantMapQualMedian  --filter-expression "MMQ < 40.0" \
            --filter-name VariantROQ           --filter-expression "ROQ < 25.0"'



rule FilterByOrientationBias:
    input:
        vcf="results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz",
        metrics="results/{token}/artifacts/{sample}/tumor_artifact.pre_adapter_detail_metrics.txt"
    output:
        "results/{token}/variantCalling/vcf/mutect2/oxog_filtered/{sample}_oxog_filtered.vcf.gz"
    params:
        reference=config["reference"],
        name="FilterByOrientationBias_{sample}",
        intervals=config["intervals_list"],
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk FilterByOrientationBias \
            -V {input.vcf} \
            --intervals {params.intervals} \
            --artifact-modes 'G/T' \
            -P {input.metrics} \
            -O {output}"
