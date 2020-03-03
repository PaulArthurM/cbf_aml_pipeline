def extra_strelka(wildcards):
    configfile: "config/config.yaml"
    extra = config['strelka']['extra']
    return extra


rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = "results/preprocessing/{sample}_G.bam",
        normal_index = "results/preprocessing/{sample}_G.bai",
        tumor = "results/preprocessing/{sample}_D.bam",
        tumor_index = "results/preprocessing/{sample}_D.bai",
        manta_candidates = "results/variantCalling/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"
    output:
        "results/variantCalling/strelka/{sample}/runWorkflow.py"
        #"results/variantCalling/strelka/{sample}/strelka_calls.vcf.gz"
    params:
        name="strelka_{sample}",
        nthread=8,
        ref = config["reference_GRCh37-lite"],
        callRegions=config['bed_intervals'],
        extra=extra_strelka
    conda:
        "../envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir results/variantCalling/strelka/{wildcards.sample} \
            --indelCandidates {input.manta_candidates} \
            --exome \
            --callRegions {params.callRegions}"



rule runWorkflow_strelka:
    input:
        "results/variantCalling/strelka/{sample}/runWorkflow.py"
    output:
        "results/variantCalling/strelka/{sample}/strelka_calls.vcf.gz"
    params:
        name="runWorkflow_strelka_{sample}",
        nthread=5,
    conda:
        "../envs/strelka.yaml"
    shell:
        "results/variantCalling/strelka/{wildcards.sample}/runWorkflow.py \
        --jobs {params.nthread} \
        -m local \
        && \
        mv results/variantCalling/strelka/{wildcards.sample}/results/variants/somatic.snvs.vcf.gz results/variantCalling/strelka/{wildcards.sample}/strelka_calls.vcf.gz"



rule mantaCandidateSmallsIndels:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam",
    output:
        #"results/variantCalling/manta/manta_{sample}.candidateSmallIndels.vcf.gz",
        "results/variantCalling/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"
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
            --runDir results/variantCalling/manta/{wildcards.sample} \
            --exome \
            --callRegions {params.callRegions} \
            && \
            results/variantCalling/manta/{wildcards.sample}/runWorkflow.py \
            --mode local \
            --jobs {params.nthread}"
