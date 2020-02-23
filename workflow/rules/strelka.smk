def extra_strelka(wildcards):
    extra = ''
    configfile: "config/config.yaml"
    if (config['strelka']['extra'] != "None"):
        for opt in config['strelka']['extra']:
            extra = extra + opt
    return extra


rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = "results/preprocessing/{sample}_G.bam",
        normal_index = "results/preprocessing/{sample}_G.bai",
        tumor = "results/preprocessing/{sample}_D.bam",
        tumor_index = "results/preprocessing/{sample}_D.bai",
        #manta_candidates = "results/variantCalling/Manta/Manta_{sample}.candidateSmallIndels.vcf.gz"
    output:
        #"results/variantCalling/strelka/{sample}/results/variants/somatic.snvs.vcf.gz"
        "results/variantCalling/strelka2/{sample}/strelka2_calls.vcf.gz"
    params:
        name="strelka_{sample}",
        nthread=8,
        ref = config["reference_GRCh37-lite"],
        extra=extra_strelka
    conda:
        "../envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir results/variantCalling/strelka2/{wildcards.sample} \
            && \
            results/variantCalling/strelka2/{wildcards.sample}/runWorkflow.py \
            --jobs {params.nthread} \
            -m local \
            {params.extra} \
            && \
            mv results/variantCalling/strelka2/{wildcards.sample}/results/variants/somatic.snvs.vcf.gz results/variantCalling/strelka2/{wildcards.sample}/strelka2_calls.vcf.gz"

"""
            --indelCandidates {input.manta_candidates} \
            --exome \
            --callRegions {params.callRegions} \

"""

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
            --runDir results/variantCalling/Manta/{wildcards.sample} \
            && \
            results/variantCalling/Manta/{wildcards.sample}/runWorkflow.py \
            --mode sge \
            --jobs {params.nthread} \
            && \
            mv results/variantCalling/Manta/{wildcards.sample}/results/variants/candidateSmallIndels.vcf.gz \
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


"""
            --exome \
            --callRegions {params.callRegions} \
"""
