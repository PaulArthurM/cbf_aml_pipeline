
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
    log:
        "logs/{token}/FilterMutectCalls/{sample}.log"
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
        --threshold-strategy OPTIMAL_F_SCORE\
        --f-score-beta 1.0 \
        -O {output}"



rule keep_pass_variants_mutect2:
    input:
        "results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.indels_snvs.vcf"
    params:
        name="keep_pass_variants_{sample}",
        nthread=5
    log:
        "logs/{token}/keep_pass_variants/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view \
        -f .,PASS \
        {input} > {output}"


# rule keep_pass_variants_strelka2_snvs:
#     input:
#         "results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.snvs.vcf.gz"
#     output:
#         "results/{token}/variantCalling/vcf/strelka2/pass/{sample}_somatic_filtered_pass.snvs.vcf"
#     params:
#         name="keep_pass_variants_{sample}",
#         nthread=5
#     log:
#         "logs/{token}/keep_pass_variants/{sample}.log"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "bcftools view \
#         -f .,PASS \
#         {input} > {output}"



# rule keep_pass_variants_strelka2_indels:
#     input:
#         "results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.indels.vcf.gz"
#     output:
#         "results/{token}/variantCalling/vcf/strelka2/pass/{sample}_somatic_filtered_pass.indels.vcf"
#     params:
#         name="keep_pass_variants_{sample}",
#         nthread=5
#     log:
#         "logs/{token}/keep_pass_variants/{sample}.log"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "bcftools view \
#         -f .,PASS \
#         {input} > {output}"


rule bcftools_merge_strelka2_vcf:
    input:
        snv="results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.snvs.vcf.gz",
        indel="results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.indels.vcf.gz"
    output:
        "results/{token}/variantCalling/vcf/{tool}/merged/{sample}.merged.vcf.gz"
    params:
        name="bcftools_merge_strelka2_vcf_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge \
        results/{wildcards.token}/variantCalling/{wildcards.tool}/{wildcards.sample}/results/variants/*vcf.gz \
        --force-samples \
        -Oz \
        -o {output}"


rule bgzip_compression:
    input:
        "results/{token}/variantCalling/vcf/{tool}/{step}/{sample}_somatic_filtered_pass.{type}.vcf"
    output:
        "results/{token}/variantCalling/vcf/{tool}/{step}/{sample}_somatic_filtered_pass.{type}.vcf.gz"
    params:
        name="bgzip_compression_{sample}",
        nthread=5
    shell:
        "bgzip {input}"


rule vcf_tabix_index:
    input:
        "results/{token}/variantCalling/vcf/{tool}/{step}/{sample}_somatic_filtered_pass.{type}.vcf.gz"
    output:
        "results/{token}/variantCalling/vcf/{tool}/{step}/{sample}_somatic_filtered_pass.{type}.vcf.gz.tbi"
    params:
        name="vcf_tabix_index_{sample}",
        nthread=5
    shell:
        "tabix {input}"


rule intersection_mutect2_strelka2_calls:
    input:
        strelka2="results/{token}/variantCalling/vcf/strelka2/pass/{sample}_somatic_filtered_pass.merged.vcf.gz",
        strelka2_tbi="results/{token}/variantCalling/vcf/strelka2/pass/{sample}_somatic_filtered_pass.merged.vcf.gz.tbi",
        mutect2="results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.indels_snvs.vcf.gz",
        mutect2_tbi="results/{token}/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.indels_snvs.vcf.gz.tbi"
    output:
        "results/{token}/variantCalling/intersection/{sample}_mutect2_strelka2_intersection.vcf"
    params:
        name="bcftools_merge_strelka2_vcf_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools isec \
        -n=2 \
        -w2 \
        {input.strelka2} {input.mutect2} > {output}"




rule bedtools_keep_diploid_variants:
    input:
        vcf="results/{token}/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz",
        diploid_regions="results/{token}/sequenza/{sample}_seqz/{sample}_segments.bed"
    output:
        diploid_snv="results/{token}/variantCalling/vcf/mutect2/diploid/{sample, [A-Za-z0-9]+}_somatic_filtered_diploid.vcf"
    params:
        name="Keep_diploid_{sample}",
        nthread=5
    log:
        "logs/{token}/bedtools_keep_diploid_variants/{sample}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.vcf} -b {input.diploid_regions} > {output}"
