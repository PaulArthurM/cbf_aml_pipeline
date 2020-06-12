rule cg_wiggle:
    input:
        ref=config['reference']
    output:
        "results/{token}/sequenza/genome_gc.wig.gz"
    params:
        name = "GC_Wiggle",
        nthread = 5,
        window = 50
    conda:
        "../envs/sequenza.yaml"
    shell:
        "sequenza-utils gc_wiggle \
             -f {input} \
             -o {output} \
             -w {params.window}"


rule sequenza_bam2seqz:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam",
        gcfile = "results/{token}/sequenza/genome_gc.wig.gz"
    output:
        gz = temp('results/{token}/sequenza/{sample}.seqz.gz'),
        tbi = temp('results/{token}/sequenza/{sample}.seqz.gz.tbi')
    params:
        reference = config['reference'],
        chrom = config['sequenza']['chrom'],
        name = "Sequenza_{sample}",
        nthread = 5
    conda:
        "../envs/sequenza.yaml"
    shell:
        "sequenza-utils bam2seqz \
            --fasta {params.reference} \
            -n {input.normal} \
            -t {input.tumor} \
            -gc {input.gcfile} \
            -o {output.gz}"


rule seqz_binning:
    input:
        'results/{token}/sequenza/{sample}.seqz.gz'
    output:
        'results/{token}/sequenza/small.{sample, [A-Za-z0-9]+}.seqz.gz'
    params:
        name="seqz_binning_{sample}",
        nthread=5
    conda:
        "../envs/sequenza.yaml"
    shell:
        "sequenza-utils seqz_binning \
        --seqz {input} \
        -w 50 \
        -o {output}"


rule sequenza_R:
    input:
        input = 'results/{token}/sequenza/small.{sample}.seqz.gz'
    output:
        output = 'results/{token}/sequenza/{sample, [A-Za-z0-9]+}_seqz/{sample}_segments.txt'
    params:
        name="Sequenza_r_{sample}",
        nthread=5,
        dir = 'results/sequenza/{sample}_seqz/'
    conda:
        "../envs/sequenza.yaml"
    script:
        "../scripts/sequenza_run.R"



rule segments_bed:
    input:
        'results/{token}/sequenza/{sample}_seqz/{sample}_segments.txt'
    output:
        'results/{token}/sequenza/{sample, [A-Za-z0-9]+}_seqz/{sample}_segments.bed'
    params:
        name = "Segments_bed_{sample}",
        nthread = 5
    shell:
        "workflow/scripts/segments_bed.sh {input} {output}"
