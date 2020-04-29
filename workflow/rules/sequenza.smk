rule cg_wiggle:
    input:
        ref=config['reference']
    output:
        "results/sequenza/genome_gc.wig.gz"
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
        gcfile = "results/sequenza/genome_gc.wig.gz"
    output:
        temp('results/sequenza/{sample}.seqz.gz')
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
            -o {output}"


rule seqz_binning:
    input:
        'results/sequenza/{sample}.seqz.gz'
    output:
        'results/sequenza/small.{sample}.seqz.gz'
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
        'results/sequenza/small.{sample}.seqz.gz'
    output:
        'results/sequenza/{sample}_seqz/Test_segments.txt'
    params:
        name="Sequenza_r_{sample}",
        nthread=5
    conda:
        "../envs/sequenza.yaml"
    script:
        "../scripts/sequenza_run.R"



rule segments_bed:
    input:
        'results/sequenza/{sample}_seqz/Test_segments.txt'
    output:
        'results/sequenza/{sample}_seqz/Test_segments.bed'
    params:
        name = "Segments_bed_{sample}",
        nthread = 5
    shell:
        "../scripts/segments_bed.sh {input} {output}"
