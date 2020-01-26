rule freebayes:
    input:
        # you can have a list of samples here
        bamTumor=config["PROJECT_DIR"] + "data/preprocessing/{sample}_D.bam",
        bamNormal=config["PROJECT_DIR"] + "data/preprocessing/{sample}_G.bam"
    output:
        config["PROJECT_DIR"] + "results/variantCalling/freebayes/raw/{sample}_freebayes.vcf"  # either .vcf or .bcf
    params:
        name="freebayes_{sample}",
        extra="",         # optional parameters
        ref=config['reference_GRCh37-lite'],
        intervals=config['intervals_list'],
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        nthreads=2
    wrapper:
        "freebayes \
            -f {params.ref} \
            --pooled-continuous \
            --pooled-discrete \
            --genotype-qualities \
            --report-genotype-likelihood-max \
            --allele-balance-priors-off \
            --min-alternate-fraction 0.03 \
            --min-repeat-entropy 1 \
            --min-alternate-count 2 \
            -t {params.intervals} \
            {input.bamTumor} \
            {input.bamNormal} > {output}"
