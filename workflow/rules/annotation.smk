rule annotate:
    input:
        "results/vcf/filtered/{sample}_somatic_filtered.vcf.gz"
    output:
        "results/vcf/annotated/{sample}_variants.funcotated.vcf"
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
        "results/variantCalling/vcf/mutect2/pass/{sample}_somatic_filtered_pass.vcf"
    output:
        "results/variantCalling/annovar/{sample}.hg19_multianno.vcf"
    params:
        path=config["PROJECT_DIR"] + "results",
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
        -out {params.path}/variantCalling/annovar/{wildcards.sample}"
