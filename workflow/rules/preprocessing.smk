# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule MarkDuplicates:
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
    shell:
        "gatk MarkDuplicates \
            -I {input} \
            -O {output.marked_bam} \
            -M {output.metrics_txt}"


#Generates recalibration table for Base Quality Score Recalibration (BQSR)
rule BaseRecalibrator:
    input:
        "results/preprocessing/{sample}_marked_duplicates.bam"
    output:
        temp("results/preprocessing/recal_data_{sample}.table")
    params:
        reference=config["reference"],
        intervals_list=config["intervals_list"],
        name="base_recalibrator_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {params.reference} \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/dbsnp_138.b37.vcf \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.b37.vcf \
            -O {output}"


#Apply base quality score recalibration
rule ApplyBQSR:
    input:
        table = "results/preprocessing/recal_data_{sample}.table",
        bam = "results/preprocessing/{sample}_marked_duplicates.bam"
    params:
        reference=config["reference"],
        name="apply_BQSR_{sample}",
        nthread=5
    output:
        temp("results/preprocessing/{sample}_marked_duplicates_BQSR.bam")
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"
