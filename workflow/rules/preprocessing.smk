# This script contain rules for:
#    - mark duplicated reads in bam files with MarkDuplicates
#    - compute a model of base-recalibration with BaseRecalibrator
#    - apply bases recalibrarion with ApplyBQSR



rule MarkDuplicates:
    """Run MarkDuplicates on un-deduplicated bam file and return a de-duplicated bam file."""
    input:
        "results/preprocessing/{sample}_{type}.{lane}.bam"
    output:
        marked_bam = temp("results/preprocessing/{sample}_{type}.{lane}_marked_duplicates.bam"),
        metrics_txt = "results/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    params:
        name="mark_duplicates_{sample}_{type}_{lane}",
        nthread=config["mark_duplicates"]["classic"]["nthread"]
    log:
        "logs/preprocessing/MarkDuplicates/{sample}_{type}.{lane}.log"
    shell:
        "gatk MarkDuplicates \
            -I {input} \
            -O {output.marked_bam} \
            -M {output.metrics_txt}"



rule BaseRecalibrator:
    """Run BaseRecalibrator on a bam file and return a model of bases-quality."""
    input:
        "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates.bam"
    output:
        temp("results/preprocessing/recal_data_{sample}_{type}.{lane}.table")
    params:
        reference=config["reference"],
        intervals_list=config["mutect2"]["intervals"],
        name="base_recalibrator_{sample}_{type}.{lane}",
        dbsnp=config["known-sites"]["dbsnp_138"],
        mills=config["known-sites"]["mills_1000G"],
        nthread=5
    log:
        "logs/preprocessing/BaseRecalibrator/{sample}_{type}.{lane}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {params.reference} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            -O {output}"



rule ApplyBQSR:
    """Run ApplyBQSR on a bam file using the associated model and create a new bam file with base's quality recalibrated."""
    input:
        table = "results/preprocessing/recal_data_{sample}_{type}.{lane}.table",
        bam = "results/preprocessing/{sample}_{type}.{lane}_marked_duplicates.bam"
    output:
        temp("results/preprocessing/{sample}_{type}.{lane, [0-9]}_marked_duplicates_BQSR.bam")
    params:
        reference=config["reference"],
        name="apply_BQSR_{sample}_{type}.{lane}",
        nthread=5
    log:
        "logs/preprocessing/ApplyBQSR/{sample}_{type}.{lane}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"
