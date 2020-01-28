def get_mutect2_input_normal_bam(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bam".format(sample = wildcards.sample)


def get_mutect2_input_tumour_bam(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bam".format(sample = wildcards.sample)


def get_mutect2_input_normal_bai(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bai".format(sample = wildcards.sample)


def get_mutect2_input_tumour_bai(wildcards):
    return config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bai".format(sample = wildcards.sample)



rule variant_calling_Mutect2:
    input:
        normal_bam="data/preprocessing/{sample}_G.bam",
        tumour_bam="data/preprocessing/{sample}_D.bam",
        normal_bai="data/preprocessing/{sample}_G.bai",
        tumour_bai="data/preprocessing/{sample}_D.bai"
    output:
        vcf_gz = config["PROJECT_DIR"] + "results/variantCalling/mutect2/raw/{sample}_mutect2.vcf.gz",
        f1r2_gz = config["PROJECT_DIR"] + "data/f1r2/{sample}_f1r2.tar.gz"
    params:
        ref=config["reference_GRCh37-lite"],
        PON=config["PON_VCF"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["intervals_list"],
        name="Mutect2_somatic_{sample}",
        nthread=config["mutect2"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk Mutect2 \
        -R {params.ref} \
        -L {params.intervals} \
        -I {input.normal_bam} \
        -I {input.tumour_bam} \
        -normal {wildcards.sample}_G_FREQEXCAP \
        --germline-resource {params.gnomad} \
        --panel-of-normals {params.PON} \
        --f1r2-tar-gz {output.f1r2_gz} \
	    --independent-mates \
        -O {output.vcf_gz}"


rule Calculate_Contamination_GetPileupSummaries:
    input:
        bam=config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bam",
	    bai=config["PROJECT_DIR"] + "data/preprocessing/{sample}_{type}.bai"
    output:
        config["PROJECT_DIR"] + "data/pileups/{sample}_{type}_pileups.table"
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
        tumour=config["PROJECT_DIR"] + "data/pileups/{sample}_D_pileups.table",
        matched=config["PROJECT_DIR"] + "data/pileups/{sample}_G_pileups.table"
    output:
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}.tumour_segmentation.tsv"
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
        config["PROJECT_DIR"] + "data/f1r2/{sample}_f1r2.tar.gz"
    output:
        config["PROJECT_DIR"] + "data/f1r2/{sample}_read-orientation-model.tar.gz"
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
        bam=config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bam",
	bai=config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bai"
    output:
        config["PROJECT_DIR"] + "data/f1r2/pileups/{sample}_D_getpileupsummaries.table"
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
        vcf=config["PROJECT_DIR"] + "results/variantCalling/mutect2/raw/{sample}_mutect2.vcf.gz",
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}.tumour_segmentation.tsv",
        orientation=config["PROJECT_DIR"] + "data/f1r2/{sample}_read-orientation-model.tar.gz"
    output:
        config["PROJECT_DIR"] + "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
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
        config["PROJECT_DIR"] + "results/variantCalling/mutect2/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        config["PROJECT_DIR"] + "results/variantCalling/mutect2/pass/{sample}_somatic_filtered_pass.vcf"
    params:
        name="keep_pass_variants_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view \
        -f .,PASS \
        {input} > {output}"
