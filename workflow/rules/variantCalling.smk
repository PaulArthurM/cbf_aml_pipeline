def extra_mutect2(wildcards):
    configfile: "config/config.yaml"
    extra = config['mutect2']['extra']
    return extra


rule variant_calling_Mutect2:
    input:
        normal_bam="results/preprocessing/{sample}_G.bam",
        tumour_bam="results/preprocessing/{sample}_D.bam",
        normal_bai="results/preprocessing/{sample}_G.bai",
        tumour_bai="results/preprocessing/{sample}_D.bai"
    output:
        vcf_gz = "results/variantCalling/mutect2/{sample}/mutect2_calls.vcf.gz",
        f1r2_gz = "results/f1r2/{sample}_f1r2.tar.gz"
    params:
        ref=config["reference_GRCh37-lite"],
        extra=extra_mutect2,
        name="Mutect2_somatic_{sample}",
        nthread=config["mutect2"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk Mutect2 \
        -R {params.ref} \
        -I {input.normal_bam} \
        -I {input.tumour_bam} \
        {params.extra} \
        -O {output.vcf_gz}"


rule Calculate_Contamination_GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
	    bai="results/preprocessing/{sample}_{type}.bai"
    output:
        "results/variantCalling/mutect2/pileups/{sample}_{type}_pileups.table"
    wildcard_constraints:
        lanes="[0-9]\.[0-9]",
    params:
        exac=config["CalculateContamination"]["GetPileupSummaries"]["exac"],
	    intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}_{type}",
        nthread=10
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.exac} \
            -L {params.intervals} \
            -O {output}"


rule Calculate_Contamination:
    input:
        tumour="results/variantCalling/mutect2/pileups/{sample}_D_pileups.table",
        matched="results/variantCalling/mutect2/pileups/{sample}_G_pileups.table"
    output:
        contamination_table="results/variantCalling/mutect2/pileups/contamination/{sample}.contamination.table",
        segmentation="results/variantCalling/mutect2/pileups/segmentation/{sample}.tumour_segmentation.tsv"
    params:
        name="CalculateContamination_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk CalculateContamination \
            -I {input.tumour} \
            -matched {input.matched} \
            --tumor-segmentation {output.segmentation} \
            -O {output.contamination_table}"

rule LearnReadOrientationModel:
    input:
        "results/variantCalling/mutect2/f1r2/{sample}_f1r2.tar.gz"
    output:
        "results/variantCalling/mutect2/f1r2/{sample}_read-orientation-model.tar.gz"
    params:
        name="LearnReadOrientationModel_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk LearnReadOrientationModel \
            -I {input} \
            -O {output}"


rule GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_D.bam",
        bai="results/preprocessing/{sample}_D.bai"
    output:
        "results/variantCalling/mutect2/f1r2/pileups/{sample}_D_getpileupsummaries.table"
    params:
        gnomad=config["mutect2"]["gnomad"]["files"]["biallelic"],
	intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}",
        nthread=5
    conda: "../envs/gatk4.yaml"
    shell:
        "gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.gnomad}  \
            -L {params.intervals} \
            -O {output}"


rule FilterMutectCalls:
    input:
        vcf="results/variantCalling/mutect2/raw/{sample}_mutect2.vcf.gz",
        contamination_table="results/variantCalling/mutect2/pileups/contamination/{sample}.contamination.table",
        segmentation="results/variantCalling/mutect2/pileups/segmentation/{sample}.tumour_segmentation.tsv",
        orientation="results/variantCalling/mutect2/f1r2/{sample}_read-orientation-model.tar.gz"
    output:
        "results/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
        #"results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    params:
        reference=config["reference_GRCh37-lite"],
        name="FilterMutectCalls_{sample}",
        nthread=config["FilterMutectCalls"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk FilterMutectCalls \
        -R {params.reference} \
        -V {input.vcf} \
        --contamination-table {input.contamination_table} \
        --tumor-segmentation {input.segmentation} \
        --orientation-bias-artifact-priors {input.orientation} \
        -O {output}"



rule keep_pass_variants:
    input:
        "results/variantCalling/vcf/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf"
        #"results/variantCalling/mutect2/pass/{sample}_somatic_filtered_pass.vcf"
    params:
        name="keep_pass_variants_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view \
        -f .,PASS \
        {input} > {output}"
