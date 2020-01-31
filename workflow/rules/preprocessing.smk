# Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
rule mark_duplicates:
    input:
        "results/preprocessing/{sample}_{type}.{lane}.bam"
    output:
        marked_bam = temp("results/preprocessing/{sample}_{type}.{lane}_marked_duplicates.bam"),
        metrics_txt = "results/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
    conda:
        "../envs/gatk4.yaml"
    params:
        name="mark_duplicates_{sample}_{type}_{lane}",
        nthread=config["mark_duplicates"]["classic"]["nthread"]
    shell:
        "gatk MarkDuplicates \
            -I {input} \
            -O {output.marked_bam} \
            -M {output.metrics_txt}"


# rule BQSRPipelineSpark:
#     input:
#         config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates.bam"
#     output:
#         temp(config["PROJECT_DIR"] + "data/bam/{sample}_marked_duplicates_BQSR.bam")
#     params:
#         reference=config["reference_GRCh37-lite"],
#         intervals_list=config["intervals_list"],
#         dbsnp_138=config["known-sites"]["dbsnp_138"],
#         mills_1000G=config["known-sites"]["mills_1000G"],
#         name="BQSRPipelineSpark_{sample}",
#         nthread=config["BQSRPipelineSpark"]["nthread"]
#     conda:
#         "../envs/gatk4.yaml"
#     shell:
#         "gatk BQSRPipelineSpark \
#         -R {params.reference} \
#         -I {input} \
#         --known-sites {params.dbsnp_138} \
#         --known-sites {params.mills_1000G} \
#         -L {params.intervals_list} \
#         -O {output} \
#         --conf 'spark.executor.cores={params.nthread}'"


#Generates recalibration table for Base Quality Score Recalibration (BQSR)
rule base_recalibrator:
    input:
        "results/preprocessing/{sample}_marked_duplicates.bam"
    output:
        temp("results/preprocessing/recal_data_{sample}.table")
    params:
        reference=config["reference_GRCh37-lite"],
        intervals_list=config["intervals_list"],
        name="base_recalibrator_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk BaseRecalibrator \
            -I {input} \
            -R {params.reference} \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/dbSNP/dbsnp_138.b37.vcf \
            --known-sites /data1/scratch/pamesl/projet_cbf/data/mills_1000G/Mills_and_1000G_gold_standard.indels.b37.vcf \
            -O {output}"


#Apply base quality score recalibration
rule apply_BQSR:
    input:
        table = "results/preprocessing/recal_data_{sample}.table",
        bam = "results/preprocessing/{sample}_marked_duplicates.bam"
    params:
        reference=config["reference_GRCh37-lite"],
        name="apply_BQSR_{sample}",
        nthread=5
    output:
        temp("results/preprocessing/{sample}_marked_duplicates_BQSR.bam")
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output}"


# # Rule for mark duplicates reads in BAM file using MarkDuplicates from GATK4
# rule mark_duplicates_spark:
#     input:
#         config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}.bam"
#     output:
#         marked_bam = temp(config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lane}_marked_duplicates.bam"),
#         metrics_txt = config["PROJECT_DIR"] + "data/metrics/{sample}_{type}.{lane}_marked_dup_metrics.txt"
#     conda:
#         "../envs/gatk4.yaml"
#     params:
#         name="mark_duplicates_spark_{sample}_{type}.{lane}",
#         nthread=config["mark_duplicates"]["spark"]["nthread"]
#     shell:
#         "gatk MarkDuplicatesSpark \
#             -I {input} \
#             -O {output.marked_bam} \
#             -M {output.metrics_txt} \
#             --conf 'spark.executor.cores={params.nthread}'"
