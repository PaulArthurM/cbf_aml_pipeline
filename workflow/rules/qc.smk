# This script contain rules for:
#    - running FastQC on bam files
#    - generating a report from all FastQC output



rule fastqc:
    """Run fastqc on input bam file and return a html report."""
    input:
        "results/preprocessing/{sample}_{type}.{lane}.bam"
    output:
        "results/{token}/quality_control/{sample}_{type}.{lane}_fastqc.html"
    params:
        name="fastq_{sample}_{type}_{lane}",
        nthread=4
    log:
        "logs/{token}/fastqc/{sample}_{type}_{lane}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input} -t {params.nthread} -o results/{wildcards.token}/quality_control/"



def getBamToQC(wildcards):
    """Return a list containing all files expected as input for MultiQC command."""
    SAMPLES = sample_sheet['samples']
    out = []
    SAMPLE = sample_sheet.set_index("samples", drop = False)
    for sample in SAMPLE["samples"]:
        for bam in SAMPLE.at[sample, "somatic_path"].split(" "):
            template = "results/{token}/quality_control/{sample}_D.{lane}_fastqc.html".format(token=wildcards.token, sample=sample, lane=getLane(bam))
            out.append(template)
        for bam in SAMPLE.at[sample, "germline_path"].split(" "):
            template = "results/{token}/quality_control/{sample}_G.{lane}_fastqc.html".format(token=wildcards.token, sample=sample, lane=getLane(bam))
            out.append(template)
    return out


rule multiqc:
    """Run multiqc on all html reports from FastQC and return a html file summarize it in one html file."""
    input:
        getBamToQC
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
        {params.input_dir} 2> {log}"
