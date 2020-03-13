rule fastqc:
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/quality_control/{sample}_{type}_fastqc.html"
    params:
        name="fastq_{sample}_{type}",
        nthread=4
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o results/quality_control/"



rule multiqc:
    input:
        expand("results/quality_control/{sample}_{type}_fastqc.html", sample=sample_sheet['samples'], type=['D', 'G'])
    output:
        "results/report/multiqc_report.html"
    params:
        name="multiqc_report",
        nthread=5,
        out_dir = "results/report/",
        input_dir = "results/quality_control/",
        out_file = "multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc \
        --force \
        -o {params.out_dir} \
        -n {params.out_file} \
        {params.input_dir}"
