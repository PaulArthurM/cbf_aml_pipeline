rule Mutect2_tumour_only:
    input:
        bam=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bam",
        bai=config["PROJECT_DIR"] + "data/bam/{sample}_G.{lane}_marked_duplicates_BQSR_merge.bai"
    output:
        temp(config["PROJECT_DIR"] + "data/vcf/{sample}_G.{lane}_marked_duplicates_BQSR_merge_for_pon.vcf.gz")
    params:
        ref=config["reference_GRCh37-lite"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["intervals_list"],
        name="Mutect2_tumour_only_{sample}_G.{lane}",
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
        "none"#VCF_IDX
    output:
        #db=directory(config["db_GDBI"]),
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
        #VCF,
        #VCF_IDX,
        test="genomicsdb.txt"
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
