# A Snakemake Pipeline for the Study of Genomic Instability

A Snakemake pipeline for the bioinformatics analysis of genomic instability from NGS data with normal/tumor pair.

</br>

<p align="center">
    <img title="Snakemake Workflow" src="https://github.com/PaulArthurM/cbf_aml_pipeline/blob/master/docs/images/Snakefile_orga.png" width=60%>
</p>

</br>

This Snakemake pipeline take aligned BAM files from germline/diagnosis samples as input and perform:

1. GATK Best Practices Preprocessing

2. Variant Calling using:

    - Mutect2 (GATK4)

    - Strelka2 + Manta

3. Ploidy/Cellularity analysis using:

    - Sequenza
    
4. (TODO) Add Williams *et al*, 2018 and/or Mobster package (Caravagna *et al*, 2019 preprint)  


