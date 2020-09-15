# An pipeline for the study of genetic instability automated with Snakemake v0.2

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
    
    - More to come.

3. Ploidy/Cellularity analysis using:

    - Sequenza
    
4. (TODO) Add Mobster package (Caravagna *et al*, 2020 Nature Genetics, Subclonal reconstruction of tumors by using machine learning and population genetics)  


