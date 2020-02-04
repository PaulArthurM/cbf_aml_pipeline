


file_sample = "/home/puissant/Sample_list.txt"
fo_sample = open(file_sample)
l = fo_sample.readlines()
samples = []
for ll in l:
    sample = ll[0:-1]
    samples.append(sample)



file = "/home/puissant/validated_SJCBF.tsv"
f_sample = open(file)
l = f_sample.readlines()
variants = []
for ll in l:
    s = ll.split("\t")
    if len(s) > 2:
        s[2] = s[2][0:-1]
        variants.append(s)


dico = {}
sign = ["KIT", "FLT3", "NRAS", "KRAS", "CBL", "JAK2"]

for sample in samples:
    dico[sample] = 0

for variant in variants:
    if variant[2] in samples:
        if variant[0] in sign:
            dico[variant[2]] += 1

print(dico)

no = 0
single = 0
int_clo = 0
for sample in dico:
    if dico[sample] == 0:
        no += 1
    if dico[sample] == 1:
        single += 1
    if dico[sample] > 1:
        int_clo += 1

for sample in dico:
    #if dico[sample] == 0:
        #print(sample)
    #if dico[sample] == 1:
        #print(sample)
    if dico[sample] > 1:
        print(sample)

if 0:
    print(no)
    print(single)
    print(int_clo)
