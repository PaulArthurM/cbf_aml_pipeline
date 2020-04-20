

rule sequenza:
    input:
        normal = "results/preprocessing/{sample}_G.bam",
        tumor = "results/preprocessing/{sample}_D.bam",
        gcfile = "/results/sequenza/genome_gc.wig.gz"
    output:
        'results/sequenza/seqzfile.{sample}.vcf'
    params:
        reference = config['reference'],
        chrom = config['sequenza']['chrom'],
        name = "Sequenza_{sample}",
        nthread = 5
    conda:
<<<<<<< HEAD
        "../envs/sequenza.yaml"
=======
        "workflow/envs/sequenza.yaml"
>>>>>>> 96699d48d257bc61c1dde95fb9ddb4562f157f44
    shell:
        "sequenza-utils.py bam2seqz \
            --fasta {params.reference} \
            -n {input.normal} \
            -t {input.tumor} \
            -gc {params.gcfile} \
            --chromosome {params.chrom} | gzip > {output}"


rule cg_wiggle:
    input:
        config['reference']
    output:
        "/results/sequenza/genome_gc.wig.gz"
    params:
        name = "GC_Wiggle",
        nthread = 5,
        window = 50
    conda:
<<<<<<< HEAD
        "../envs/sequenza.yaml"
=======
        "workflow/envs/sequenza.yaml"
>>>>>>> 96699d48d257bc61c1dde95fb9ddb4562f157f44
    shell:
        "sequenza-utils gc_wiggle \
            -f {input} \
            -O {output} \
            -w {params.window}"
