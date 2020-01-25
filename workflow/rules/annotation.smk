rule annotate:
    input:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered.vcf.gz"
    output:
        config["PROJECT_DIR"] + "data/vcf/annotated/{sample}_{normal_lanes}-{tumour_lanes}_variants.funcotated.vcf"
    params:
        reference=config["reference_GRCh37-lite"],
        data_sources=config["funcotator"]["directory"],
        name="Funcotator_annotate_{sample}",
        nthread=4
    conda:
        "../envs/gatk4.yaml"
    shell:
        "gatk Funcotator \
        --variant {input} \
        --reference {params.reference} \
        --ref-version b37 \
        --data-sources-path {params.data_sources} \
        --output {output} \
        --output-file-format VCF"



rule annovar:
    input:
        config["PROJECT_DIR"] + "data/vcf/filtered/{sample}_{normal_lanes}-{tumour_lanes}_somatic_filtered_pass.vcf"
    output:
        config["PROJECT_DIR"] + "data/vcf/annotated/{sample}_{normal_lanes}-{tumour_lanes}.avinput"
    params:
        path="/data1/scratch/pamesl/projet_cbf/data",
        name="Annovar_annotate_{sample}",
        nthread=4
    shell:
        "perl /data1/scratch/pamesl/projet_cbf/annovar/table_annovar.pl \
        {input} \
        /data1/scratch/pamesl/projet_cbf/annovar/humandb \
        -protocol refGene \
        -operation g \
        -nastring '.' \
        --buildver hg19 \
        --remove \
        --vcfinput \
        -out {params.path}/vcf/annotated/{wildcards.sample}_{wildcards.normal_lanes}-{wildcards.tumour_lanes}"
