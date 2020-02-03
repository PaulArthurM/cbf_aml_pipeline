rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = "results/preprocessing/{sample}_G.bam",
        normal_index = "results/preprocessing/{sample}_G.bai",
        tumor = "results/preprocessing/{sample}_D.bam",
        tumor_index = "results/preprocessing/{sample}_D.bai",
        manta_candidates = "results/variantCalling/Manta/Manta_{sample}.candidateSmallIndels.vcf.gz"
    output:
        "results/variantCalling/Strelka/Strelka_{sample}_variants.vcf.gz"
    params:
        name="strelka_{sample}",
        nthread=8,
        ref = config["reference_GRCh37-lite"],
        callRegions = config["bed_intervals"]
    conda:
        "../envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --indelCandidates {input.manta_candidates} \
            --exome \
            --callRegions {params.callRegions} \
            --runDir results/variantCalling/Strelka/{wildcards.sample} \
            && \
            results/variantCalling/Strelka/{wildcards.sample}/runWorkflow.py \
            --mode sge \
            --jobs {params.nthread} \
            && \
            mv results/variantCalling/Strelka/{wildcards.sample}/results/variants/genome.*.vcf.gz \
                results/variantCalling/Strelka/Strelka_{wildcards.sample}_genome.vcf.gz \
            mv results/variantCalling/Strelka/{wildcards.sample}/results/variants/genome.*.vcf.gz.tbi \
                results/variantCalling/Strelka/Strelka_{wildcards.sample}_genome.vcf.gz.tbi \
            mv results/variantCalling/Strelka/{wildcards.sample}/results/variants/variants.vcf.gz \
                results/variantCalling/Strelka/Strelka_{wildcards.sample}_variants.vcf.gz \
            mv results/variantCalling/Strelka/{wildcards.sample}/results/variants/variants.vcf.gz.tbi \
                results/variantCalling/Strelka/Strelka_{wildcards.sample}_variants.vcf.gz.tbi"


rule mantaCandidateSmallsIndels:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam",
    output:
        "results/variantCalling/Manta/Manta_{sample}.candidateSmallIndels.vcf.gz",
    params:
        name="Manta_{sample}",
        nthread=8,
        ref = config["reference_GRCh37-lite"],
        callRegions = config["bed_intervals"]
    conda:
        "../envs/strelka.yaml"
    shell:
        "configManta.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --exome \
            --callRegions {params.callRegions} \
            --runDir results/variantCalling/Manta/{wildcards.sample} \
            && \
            results/variantCalling/Manta/{wildcards.sample}/runWorkflow.py \
            --mode sge \
            --jobs {params.nthread} \
            && \
<<<<<<< HEAD
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/candidateSmallIndels.vcf.gz \
=======
            mv results/variantCalling/Manta/{wildcards.ample}/results/variants/candidateSmallIndels.vcf.gz \
>>>>>>> df7da8f685331be6f4596f46e7f734cad0590773
                results/variantCalling/Manta/Manta_{wildcards.sample}.candidateSmallIndels.vcf.gz \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/candidateSmallIndels.vcf.gz.tbi \
                results/variantCalling/Manta/Manta_{wildcards.sample}.candidateSmallIndels.vcf.gz.tbi \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/candidateSV.vcf.gz \
                results/variantCalling/Manta/Manta_{wildcards.sample}.candidateSV.vcf.gz \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/candidateSV.vcf.gz.tbi \
                results/variantCalling/Manta/Manta_{wildcards.sample}.candidateSV.vcf.gz.tbi \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/tumorSV.vcf.gz \
                results/variantCalling/Manta/Manta_{wildcards.sample}.tumorSV.vcf.gz \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/tumorSV.vcf.gz.tbi \
                results/variantCalling/Manta/Manta_{wildcards.sample}.tumorSV.vcf.gz.tbi"
