# This script contain rules for:
#    - running FastQC on bam files
#    - generating a report from all FastQC output



rule fastqc:
    """Run fastqc on input bam file and return a html report."""
    input:
        "results/preprocessing/{sample}_{type}.bam"
    output:
        "results/{token}/quality_control/{sample}_{type}_fastqc.html"
    params:
        name="fastq_{sample}_{type}",
        nthread=4
    log:
        "logs/{token}/fastqc/{sample}_{type}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o results/{wildcards.token}/quality_control/"



rule multiqc:
    """Run multiqc on all html reports from FastQC and return a html file summarize it in one html file."""
    input:
        expand("results/{token}/quality_control/{sample}_{type}_fastqc.html", sample=sample_sheet['samples'], type=['D', 'G'], token = config['token'])
    output:
        "results/{token}/quality_control/report/multiqc_report.html"
    params:
        name="multiqc_report",
        nthread=5,
        out_dir = "results/{token}/quality_control/report/",
        input_dir = "results/{token}/quality_control/",
        out_file = "multiqc_report.html"
    log:
        "logs/{token}/multiqc/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc \
        --force \
        -o {params.out_dir} \
        -n {params.out_file} \
        {params.input_dir}"
