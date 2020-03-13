from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

project = "SJCBF-AML"

genome = "GRCh37"

vcfFiles = "/home/puissant/old_cbf_aml_pipeline/data/vcf/mutect2"

matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, vcfFiles, exome=True, bed_file=None, chrom_based=False, plot=True, tsb_stat=True, seqInfo=True)
