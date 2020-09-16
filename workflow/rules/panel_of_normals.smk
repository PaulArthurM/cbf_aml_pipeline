rule Mutect2_tumour_only:
    input:
        bam="results/preprocessing/{sample}_G.bam",
        bai="results/preprocessing/{sample}_G.bai"
    output:
        "results/{token}/vcf/{sample}_pon.vcf.gz"
    params:
        ref=config["reference"],
        gnomad=config["mutect2"]["gnomad"]["files"]["raw"],
        intervals=config["mutect2"]["intervals"],
        name="Mutect2_tumour_only_{sample}",
        nthread=config["mutect2"]["nthread"]
    log:
        "logs/{token}/Mutect2_tumour_only/{sample}.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        " gatk Mutect2 \
        -R {params.ref} \
        -I {input.bam} \
        -max-mnp-distance 0 \
        -L {params.intervals} \
        -O {output} 2> {log}"


rule GenomicsDBImport:
    input:
        expand("results/{token}/vcf/{sample}_pon.vcf.gz", token=config["token"], sample=pep.sample_table["sample_name"])
    output:
        db=directory("results/{token}/GenomicsDBImport".format(token=config["token"])),
        test="results/{token}/genomicsdb.txt".format(token=config["token"])
    params:
        ref=config["reference"],
        inputString = lambda wildcards, input: " -V ".join(input),
        intervals=config["mutect2"]["intervals"],
        name="GenomicsDB",
        nthread=20
    log:
        "logs/{token}/GenomicsDBImport/logs.log".format(token=config["token"])
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk GenomicsDBImport \
        -R {params.ref} \
        -L {params.intervals} \
        --genomicsdb-workspace-path {output.db} \
        --merge-input-intervals \
        -V {params.inputString} && touch {output.test}"


rule CreateSomaticPanelOfNormals:
    input:
        "results/{token}/genomicsdb.txt"
    output:
        "results/{token}/pon/panel_of_normals.vcf.gz"
    params:
        ref=config["reference"],
        db= config["directory"] + "results/{token}/GenomicsDBImport",
        name="create_PON",
        nthread=20
    log:
        "logs/{token}/CreateSomaticPanelOfNormals/logs.log"
    conda:
        "../envs/gatk4.1.7.0.yaml"
    shell:
        "gatk CreateSomaticPanelOfNormals \
        -R {params.ref} \
        -V gendb://{params.db} \
        -O {output} 2> {log}"
