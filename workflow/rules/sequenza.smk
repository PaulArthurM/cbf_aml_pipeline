


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
        "../envs/sequenza.yaml"
    shell:
        "sequenza-utils.py bam2seqz \
            --fasta {params.reference} \
            -n {input.normal} \
            -t {input.tumor} \
            -gc {params.gcfile} \
            --chromosome {params.chrom} | gzip > {output}"


<<<<<<< HEAD
rule cg_wiggle:
    input:
        ref=config['reference']
    output:
        "/results/sequenza/genome_gc.wig.gz"
    params:
        name = "GC_Wiggle",
        nthread = 5,
        window = 50
    conda:
        "../envs/sequenza.yaml"
    shell:
        "sequenza-utils gc_wiggle \
=======
 rule cg_wiggle:
     input:
         ref=config['reference']
     output:
         "/results/sequenza/genome_gc.wig.gz"
     params:
         name = "GC_Wiggle",
         nthread = 5,
         window = 50
     conda:
         "../envs/sequenza.yaml"
     shell:
         "sequenza-utils gc_wiggle \
>>>>>>> 11c3a94f58ee67c658e3f58403917439dff017d3
             -f {input} \
             -O {output} \
             -w {params.window}"
