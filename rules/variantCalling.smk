rule variant_calling_Mutect2:
    input:
        normal_bam = config["PROJECT_DIR"] + "data/bam/{sample}_G.{lanes_normal}_marked_duplicates_BQSR_merge.bam",
        tumour_bam = config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes_tumour}_marked_duplicates_BQSR_merge.bam",
	    normal_bai = config["PROJECT_DIR"] + "data/bam/{sample}_G.{lanes_normal}_marked_duplicates_BQSR_merge.bai",
	    tumour_bai = config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes_tumour}_marked_duplicates_BQSR_merge.bai"
    output:
        vcf_gz = config["PROJECT_DIR"] + "data/vcf/{sample}_{lanes_normal}-{lanes_tumour}_somatic.vcf.gz",
        f1r2_gz = config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_f1r2.tar.gz"
    wildcard_constraints:
        lanes_tumour="[0-9]\.[0-9]",
        lanes_normal="[0-9]\.[0-9]",
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
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lanes}_marked_duplicates_BQSR_merge.bam",
	    bai=config["PROJECT_DIR"] + "data/bam/{sample}_{type}.{lanes}_marked_duplicates_BQSR_merge.bai"
    output:
        config["PROJECT_DIR"] + "data/pileups/{sample}_{type}.{lanes}_pileups.table"
    wildcard_constraints:
        lanes="[0-9]\.[0-9]",
    params:
        exac=config["CalculateContamination"]["GetPileupSummaries"]["exac"],
	intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}_{type}.{lanes}",
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
        tumour=config["PROJECT_DIR"] + "data/pileups/{sample}_D.{tumour_lanes}_pileups.table",
        matched=config["PROJECT_DIR"] + "data/pileups/{sample}_G.{normal_lanes}_pileups.table"
    output:
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}_{normal_lanes}-{tumour_lanes}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}_{normal_lanes}-{tumour_lanes}.tumour_segmentation.tsv"
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
        config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_f1r2.tar.gz"
    output:
        config["PROJECT_DIR"] + "data/f1r2/{sample}_{lanes_normal}-{lanes_tumour}_read-orientation-model.tar.gz"
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
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes}_marked_duplicates_BQSR_merge.bam",
	bai=config["PROJECT_DIR"] + "data/bam/{sample}_D.{lanes}_marked_duplicates_BQSR_merge.bai"
    output:
        config["PROJECT_DIR"] + "data/f1r2/pileups/{sample}_D.{lanes}_getpileupsummaries.table"
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
        vcf=config["PROJECT_DIR"] + "data/vcf/{sample}_{normal_lanes}-{tumour_lanes}_somatic.vcf.gz",
        contamination_table=config["PROJECT_DIR"] + "data/pileups/contamination/{sample}_{normal_lanes}-{tumour_lanes}.contamination.table",
        segmentation=config["PROJECT_DIR"] + "data/pileups/segmentation/{sample}_{normal_lanes}-{tumour_lanes}.tumour_segmentation.tsv",
        orientation=config["PROJECT_DIR"] + "data/f1r2/{sample}_{normal_lanes}-{tumour_lanes}_read-orientation-model.tar.gz"
    output:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered.vcf.gz"
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
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered.vcf.gz"
    output:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered_pass.vcf.gz"
    params:
        name="keep_pass_variants_{sample}",
        nthread=5
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view \
        -f .,PASS \
        {input} > {output}"
