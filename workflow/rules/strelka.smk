rule strelka:
    input:
        # The normal bam and its index
        # are optional input
        normal = "results/preprocessing/{sample}_G.bam",
        normal_index = "results/preprocessing/{sample}_G.bai",
        tumor = "results/preprocessing/{sample}_D.bam",
        tumor_index = "results/preprocessing/{sample}_D.bai",
        manta_candidates = "results/{token}/variantCalling/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"
    output:
        "results/{token}/variantCalling/strelka/{sample}/runWorkflow.py"
    params:
        name="strelka_{sample}",
        nthread=8,
        ref = config["reference"],
        callRegions=config['bed_intervals']
        #extra=extra_strelka
    log:
        "logs/{token}/strelka/{sample}.log"
    conda:
        "../envs/strelka.yaml"
    shell:
        "configureStrelkaSomaticWorkflow.py  \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir results/{wildcards.token}/variantCalling/strelka/{wildcards.sample} \
            --indelCandidates {input.manta_candidates} \
            --exome \
            --callRegions {params.callRegions}"



rule runWorkflow_strelka:
    input:
        "results/{token}/variantCalling/strelka/{sample}/runWorkflow.py"
    output:
        snv="results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.snvs.vcf",
        indel="results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.indels.vcf"
    params:
        name="runWorkflow_strelka_{sample}",
        nthread=5
    log:
        "logs/{token}/runWorkflow_strelka/{sample}.log"
    conda:
        "../envs/strelka.yaml"
    shell:
        "results/{wildcards.token}/variantCalling/strelka/{wildcards.sample}/runWorkflow.py \
        --jobs {params.nthread} \
        -m local && \
        mv  results/{wildcards.token}/variantCalling/strelka/{wildcards.sample}/results/variants/somatic.snvs.vcf.gz results/{wildcards.token}/variantCalling/strelka/{wildcards.sample}/results/variants/somatic.snvs.vcf && \
        mv results/{wildcards.token}/variantCalling/strelka/{wildcards.sample}/results/variants/somatic.indels.vcf.gz results/{wildcards.token}/variantCalling/strelka/{wildcards.sample}/results/variants/somatic.indels.vcf"


rule strelka2_bgzip_compression:
    input:
        "results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.{type}.vcf"
    output:
        "results/{token}/variantCalling/strelka/{sample}/results/variants/somatic.{type}.vcf.gz"
    params:
        name="strelka2_bgzip_compression_{sample}",
        nthread=5
    shell:
        "bgzip {input}"


rule mantaCandidateSmallsIndels:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam",
    output:
        "results/{token}/variantCalling/manta/{sample}/results/variants/candidateSmallIndels.vcf.gz"
    params:
        name="Manta_{sample}",
        nthread=8,
        ref = config["reference"],
        callRegions = config["bed_intervals"]
    log:
        "logs/{token}/mantaCandidateSmallsIndels/{sample}.log"
    conda:
        "../envs/strelka.yaml"
    shell:
        "configManta.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.ref} \
            --runDir results/{wildcards.token}/variantCalling/manta/{wildcards.sample} \
            --exome \
            --callRegions {params.callRegions} \
            && \
            results/{wildcards.token}/variantCalling/manta/{wildcards.sample}/runWorkflow.py \
            --mode local \
            --jobs {params.nthread}"
