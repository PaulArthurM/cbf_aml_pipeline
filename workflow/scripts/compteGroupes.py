


file_sample = "/data1/scratch/pamesl/cbf_aml_pipeline/Sample_File_SJCBF.map"
fo_sample = open(file_sample)
samples = fo_sample.readlines()

file = "/data1/scratch/pamesl/cbf_aml_pipeline/validated_SJCBF.tsv"
f_sample = open(file)
variants = f_sample.readlines()


dico = {}

for variant in variants:
    sample = variant.split("\t")[2]
    if sample in samples:
        if sample in dico:
            dico[sample] += 1
        elif not sample in dico:
            dico[sample] = 1


print(dico)
