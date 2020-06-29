# Yet Another Pipeline for Normal/Tumor Variant Calling from Whole-Exome Sequencing data v0.2

A Snakemake pipeline for the bioinformatics analysis of WES data with normal/tumor pair.

</br>

<p align="center">
    <img title="Snakemake Workflow" src="https://github.com/PaulArthurM/cbf_aml_pipeline/blob/master/docs/images/Snakefile_orga.png" width=60%>
</p>

</br>

This Snakemake's pipeline take aligned BAM files from germline/diagnosis samples as input and perform:

1. GATK Best Practices Preprocessing

2. Variant Calling using:

    - Mutect2 (GATK4)

    - Strelka2 + Manta

3. Ploidy/Cellularity analysis using:

    - Sequenza
    
4. (TODO) Add Williams *et al*, 2018 and/or Mobster package (Caravagna *et al*, 2019 preprint)  


