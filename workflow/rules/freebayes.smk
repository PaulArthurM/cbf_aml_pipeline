rule freebayes:
    input:
        # you can have a list of samples here
        bamTumor="results/preprocessing/{sample}_D.bam",
        bamNormal="results/preprocessing/{sample}_G.bam"
    output:
        "results/{token}/variantCalling/freebayes/{sample}/freebayes_calls.vcf"  # either .vcf or .bcf
    params:
        name="freebayes_{sample}",
        #extra=extra_freebayes,         # optional parameters
        ref=config['reference'],
        intervals=config['bed_intervals'],
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        nthread=5
    log:
        "logs/{token}/freebayes/{sample}.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes \
            -f {params.ref} \
            -r X \
            --pooled-continuous \
            --pooled-discrete \
            --genotype-qualities \
            --report-genotype-likelihood-max \
            --allele-balance-priors-off \
            --min-alternate-fraction 0.05 \
            --min-repeat-entropy 1 \
            --min-alternate-count 2 \
            {input.bamTumor} \
            {input.bamNormal} > {output}"
