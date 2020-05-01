"""
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

project = "SJCBF-AML"

genome = "GRCh37"

vcfFiles = "/home/puissant/old_cbf_aml_pipeline/data/vcf/mutect2"

matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, vcfFiles, exome=True, bed_file=None, chrom_based=False, plot=True, tsb_stat=True, seqInfo=True)
"""

sample = snakemake.wildcards["sample"]
genome_version = snakemake.params["genome_version"]
out_dir = snakemake.params["out_dir"]
bed = snakemake.input["bed"]


from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37', rsync=False, bash=True)


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc(sample, genome_version, out_dir, plot=True, exome=True, bed_file=bed, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


matrix_path = out_dir + "output/SBS/{sample}.SBS96.exome".format(sample=sample)

from sigProfilerPlotting import sample_portrait as sP
sP.samplePortrait(sample_matrices_path, output_path, project, percentage=False)
