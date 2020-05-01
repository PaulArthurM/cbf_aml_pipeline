rule variant_calling_Mutect2:
    input:
        normal_bam="results/preprocessing/{sample}_G.bam",
        tumour_bam="results/preprocessing/{sample}_D.bam",
        normal_bai="results/preprocessing/{sample}_G.bai",
        tumour_bai="results/preprocessing/{sample}_D.bai"
    output:
        vcf_gz = "results/{token}/variantCalling/mutect2/{sample}/mutect2_calls.vcf.gz",
        f1r2_gz = "results/{token}/variantCalling/mutect2/f1r2/{sample}_f1r2.tar.gz"
    params:
        ref=config["reference"],
        extra=extra_mutect2,
        name="Mutect2_somatic_{sample}",
        nthread=config["mutect2"]["nthread"]
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
        {params.extra} \
        -O {output.vcf_gz}"


rule Calculate_Contamination_GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_{type}.bam",
	    bai="results/preprocessing/{sample}_{type}.bai"
    output:
        "results/{token}/variantCalling/mutect2/pileups/{sample}_{type}_pileups.table"
    wildcard_constraints:
        lanes="[0-9]\.[0-9]",
    params:
        exac=config["CalculateContamination"]["GetPileupSummaries"]["exac"],
	    intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}_{type}",
        nthread=10
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        " gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.exac} \
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
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk LearnReadOrientationModel \
            -I {input} \
            -O {output}"


rule GetPileupSummaries:
    input:
        bam="results/preprocessing/{sample}_D.bam",
        bai="results/preprocessing/{sample}_D.bai"
    output:
        "results/{token}/variantCalling/mutect2/f1r2/pileups/{sample}_D_getpileupsummaries.table"
    params:
        gnomad=config["mutect2"]["gnomad"]["files"]["biallelic"],
	    intervals=config["intervals_list"],
        name="GetPileupSummaries_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.gnomad}  \
            -L {params.intervals} \
            -O {output}"



rule CollectSequencingArtifactMetrics:
    input:
        "results/preprocessing/{sample}_D.bam"
    output:
        "results/{token}/artifacts/{sample}/tumor_artifact.pre_adapter_detail_metrics.txt"
    params:
        reference=config["reference"],
        path_out="results/{token}/artifacts/{sample}/tumor_artifact",
        name="CollectSequencingArtifactMetrics_{sample}",
        nthread=5
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        ' gatk CollectSequencingArtifactMetrics \
            -R {params.reference} \
            -I {input} \
            --FILE_EXTENSION ".txt" \
            -O {params.path_out}'
