rule Mutect2_tumour_only:
    input:
        bam="results/preprocessing/{sample}_G.bam",
        bai="results/preprocessing/{sample}_G.bai"
    output:
        temp("results/vcf/{sample}_pon.vcf.gz")
    params:
        ref=config["reference_GRCh37-lite"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["intervals_list"],
        name="Mutect2_tumour_only_{sample}",
        nthread=config["mutect2"]["nthread"]
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk Mutect2 \
        -R {params.ref} \
        -I {input.bam} \
        -max-mnp-distance 0 \
        -L {params.intervals} \
        -O {output}"


rule GenomicsDB:
    input:
        PON_VCF
    output:
        db=directory(config["db_GDBI"]),
        test="genomicsdb.txt"
    params:
        ref=config["reference_GRCh37-lite"],
        inputString = lambda wildcards, input: " -V ".join(input),
        intervals=config["intervals_list"],
        name="GenomicsDB",
        nthread=20
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk GenomicsDBImport \
        -R {params.ref} \
        -L {params.intervals} \
        --genomicsdb-workspace-path {output.db} \
        --merge-input-intervals \
        -V {params.inputString} && touch {output.test}"


rule CreateSomaticPanelOfNormals:
    input:
        "genomicsdb.txt"
    output:
        config["PON_VCF"]
    params:
        ref=config["reference_GRCh37-lite"],
        db=config["db_GDBI"],
        name="create_PON",
        nthread=20
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk CreateSomaticPanelOfNormals \
        -R {params.ref} \
        -V gendb://{params.db} \
        -O {output}"
