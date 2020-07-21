# This script contain rules for:
#    - call candidates variants with Mutect2
#    - compute contamination for each sample
#    - learn orientation model for artifacts filtering



rule variant_calling_Mutect2:
    input:
        normal_bam="results/preprocessing/{sample}_G.bam",
        tumour_bam="results/preprocessing/{sample}_D.bam",
        normal_bai="results/preprocessing/{sample}_G.bai",
        tumour_bai="results/preprocessing/{sample}_D.bai",
        #pon="results/{token}/pon/panel_of_normals.vcf.gz"
    output:
        vcf_gz = "results/{token}/variantCalling/mutect2/{sample}/mutect2_calls.vcf.gz",
        f1r2_gz = "results/{token}/variantCalling/mutect2/f1r2/{sample}_f1r2.tar.gz"
    params:
        ref=config["reference"],
        gnomad_raw=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["mutect2"]["intervals"],
        name="Mutect2_somatic_{sample}",
        nthread=config["mutect2"]["nthread"]
    log:
        "logs/{token}/variant_calling_Mutect2/{sample}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk Mutect2 \
        -R {params.ref} \
        -I {input.normal_bam} \
        -I {input.tumour_bam} \
        --tumor {wildcards.sample}_D_FREQEXCAP \
        --normal {wildcards.sample}_G_FREQEXCAP \
        --f1r2-tar-gz {output.f1r2_gz} \
        --germline-resource {params.gnomad_raw} \
        --L {params.intervals} \
        -O {output.vcf_gz}"

#        --panel-of-normals {input.pon} \


rule Calculate_Contamination_GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
	    bai="results/preprocessing/{sample}_{type}.bai",
        exac="ressources/GetPileupSummaries/somatic-b37_small_exac_common_3.vcf"
    output:
        "results/{token}/variantCalling/mutect2/pileups/{sample}_{type}_pileups.table"
    params:
        intervals=config["mutect2"]["intervals"],
        name="GetPileupSummaries_{sample}_{type}",
        nthread=10
    log:
        "logs/{token}/Calculate_Contamination_GetPileupSummaries/{sample}_{type}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        " gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.exac} \
            -L {params.intervals} \
            -O {output}"


rule Calculate_Contamination:
    input:
        tumour="results/{token}/variantCalling/mutect2/pileups/{sample}_D_pileups.table",
        matched="results/{token}/variantCalling/mutect2/pileups/{sample}_G_pileups.table"
    output:
        contamination_table="results/{token}/variantCalling/mutect2/pileups/contamination/{sample}.contamination.table",
        segmentation="results/{token}/variantCalling/mutect2/pileups/segmentation/{sample}.tumour_segmentation.tsv"
    params:
        name="CalculateContamination_{sample}",
        nthread=5
    log:
        "logs/{token}/CalculateContamination/{sample}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        " gatk CalculateContamination \
            -I {input.tumour} \
            -matched {input.matched} \
            --tumor-segmentation {output.segmentation} \
            -O {output.contamination_table}"


rule LearnReadOrientationModel:
    input:
        "results/{token}/variantCalling/mutect2/f1r2/{sample}_f1r2.tar.gz"
    output:
        "results/{token}/variantCalling/mutect2/f1r2/{sample}_read-orientation-model.tar.gz"
    params:
        name="LearnReadOrientationModel_{sample}",
        nthread=5
    log:
        "logs/{token}/LearnReadOrientationModel/{sample}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk LearnReadOrientationModel \
            -I {input} \
            -O {output}"


rule GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_D.bam",
        bai="results/preprocessing/{sample}_D.bai",
        gnomad_biallelic="ressources/somatic-b37_af-only-gnomad.b37.BIALLELIC.vcf"
    output:
        "results/{token}/variantCalling/mutect2/f1r2/pileups/{sample}_D_getpileupsummaries.table"
    params:
        name="GetPileupSummaries_{sample}",
        nthread=5
    log:
        "logs/{token}/GetPileupSummaries/{sample}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad_biallelic}  \
            -L {input.gnomad_biallelic} \
            -O {output}"
