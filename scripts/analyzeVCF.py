import re
import argparse


class Sample():
    def __init__(self, sample_name, variants):

        CLONAL_INTERFERENCE_GENES = ["KIT", "FLT3", "NRAS", "KRAS", "JAK2", "CBL"]  # temporaire


        def getClonalInterferenceVariants(variants):
            cpt  = 0
            genes = []
            nTotal = 0
            for variant in variants:
                nTotal += 1
                if variant.geneName in CLONAL_INTERFERENCE_GENES:
                    cpt+=1
                    genes.append(variant.geneName)
            return [cpt, genes, nTotal]


        def getCategorie(nClonalInterferenceVariants):
            if nClonalInterferenceVariants == 0:
                return "Group No Mutation"
            elif nClonalInterferenceVariants == 1:
                return "Group Single Mutation"
            elif nClonalInterferenceVariants >= 2:
                return "Group Clonal Interference"
            else:
                return "ERROR: No group."


        self.sample_name = sample_name
        self.nClonalInterferenceVariants, self.clonalInterferenceVariants, self.nVariantTotal = getClonalInterferenceVariants(variants)
        self.categorie = getCategorie(self.nClonalInterferenceVariants)




class Variant():
    def __init__(self, line, vcf_file):

        def get_sample(vcf_file):
            m = re.search("(SJCBF[0-9]+)", vcf_file)
            if m:
                return m.group(1)


        def get_chr(line):
            return line.split("\t")[0]


        def get_pos(line):
            return line.split("\t")[1]


        def get_ref(line):
            return line.split("\t")[3]


        def get_alt(line):
            return line.split("\t")[4]


        def isExonic(line):
            m = re.search("Func.refGene=exonic", line)
            if m:
                return True
            else:
                return False


        def get_geneName(line):
            m = re.search("Gene.refGene=([A-Za-z0-9]+);", line)
            if m:
                return m.group(1)

        def get_exonicFunc(line):
            m = re.search("ExonicFunc.refGene=([a-zA-z0-9]+);", line)
            if m:
                return m.group(1)


        self.sample = get_sample(vcf_file)
        self.chr = get_chr(line)
        self.pos = get_pos(line)
        # if isExonic(line):
        self.ref = get_ref(line)
        self.alt = get_alt(line)
        self.geneName = get_geneName(line)
        self.exonicFunc = get_exonicFunc(line)
        # else:
        #     self.ref = "Not exonic"
        #     self.alt = "Not exonic"
        #     self.geneName = "Not exonic"
        #     self.exonicFunc = "Not exonic"


def get_sample(vcf_file):
    m = re.search("(SJCBF[0-9]+)", vcf_file)
    if m:
        return m.group(1)

def open_file(file_path):
    if file_path[-1] == "\n":
        file = open(file_path[0:-1], 'r')
    else:
        file = open(file_path, 'r')
    return file.readlines()


def readVCF(vcf_file):
    lines = open_file(vcf_file)
    return [get_sample(vcf_file), [Variant(line, vcf_file) for line in lines if line[0] != "#"]]


def showAllSamplesInfo(samples):
    print('#SAMPLE\tCHROM\tPOS\tREF\tALT\tGENE_NAME\tFUNC\n')
    for sample in samples:
        for variant in samples[sample]:
            txt = "{sample_name}\t{chr}\t{pos}\t{ref}\t{alt}\t{geneName}\t{exonicFunc}\n".format(sample_name=sample, geneName=variant.geneName, exonicFunc=variant.exonicFunc, chr=variant.chr, pos=variant.pos, ref=variant.ref, alt=variant.alt)
            print(txt)


def assignCategorie(samples):
    grouped_samples = {'Group No Mutation':[], 'Group Single Mutation':[], 'Group Clonal Interference':[]}
    for sample in samples:
        instance = Sample(sample, samples[sample])
        grouped_samples[instance.categorie].append(instance)
    return grouped_samples


def showGroupedSamples(grouped_samples):
    print("#GROUP\tN_CI\tN_TOT\tSAMPLE")
    for group in grouped_samples:
        for sample in grouped_samples[group]:
            txt = "{group}\t{nClonalInterferenceVariants}\t{nTot}\t{sample}"
            print(txt.format(nTot=sample.nVariantTotal, sample=sample.sample_name, nClonalInterferenceVariants=sample.nClonalInterferenceVariants, group=sample.categorie))



def main(args):
    vcfs = open_file(args.v)
    samples = {}
    for vcf in vcfs:
        k,v = readVCF(vcf)
        samples[k] = v
    showAllSamplesInfo(samples)
    #grouped_samples = assignCategorie(samples)
    #showGroupedSamples(grouped_samples)


if __name__== '__main__':
    parser = argparse.ArgumentParser(description='Analyze VCF.')
    parser.add_argument('-v', default=None, required=True, type=str, help="File with all VCFs name.")
    # parser.add_argument('-e', default='SJCBF', type=str, help="Experience. Default: SJCBF.")
    # parser.add_argument('-j', default=True, type=bool, help="Create JSON. Default: True")
    # parser.add_argument('-p', default='/data1/scratch/pamesl/projet_cbf/data/bam/', type=str, help="Path to bam files repository.")
    # parser.add_argument('-l', default=1, type=int, help="Number of files to download. Default: 1.")
    # parser.add_argument('-t', default=None, type=str, help="Type of files to download. For exemple: 'D' or 'G'.", required=True)

    args = parser.parse_args()

    main(args)
