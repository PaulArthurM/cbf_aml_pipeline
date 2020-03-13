rule freebayes:
    input:
        # you can have a list of samples here
        bamTumor="results/preprocessing/{sample}_D.bam",
        bamNormal="results/preprocessing/{sample}_G.bam"
    output:
        "results/variantCalling/freebayes/{sample}/freebayes_calls.vcf"  # either .vcf or .bcf
    params:
        name="freebayes_{sample}",
        extra=extra_freebayes,         # optional parameters
        ref=config['reference'],
        intervals=config['bed_intervals'],
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        nthread=5
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes \
            -f {params.ref} \
            -t {params.intervals} \
            {input.bamTumor} \
            {input.bamNormal} > {output}"
