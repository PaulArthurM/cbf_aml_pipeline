import re
import argparse


class Variant():
    def __init__(self, line, vcf_file):

        def get_sample(vcf_file):
            m = re.search("(SJCBF[0-9]+)", vcf_file)
            if m:
                return m.group(1)

        def get_ref(line):
            return line.split("\t")[2]

        def get_alt(line):
            return line.split("\t")[3]

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
        if isExonic(line):
            self.ref = get_ref(line)
            self.alt = get_alt(line)
            self.geneName = get_geneName(line)
            self.exonicFunc = get_exonicFunc(line)
        else:
            self.ref = "Not exonic"
            self.alt = "Not exonic"
            self.geneName = "Not exonic"
            self.exonicFunc = "Not exonic"


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
    return {get_sample(vcf_file):[Variant(line, vcf_file) for line in lines if line[0] != "#"]}


def showAllSamplesInfo(samples):
    for sample in samples.keys():
        print(sample)
    # for variant in variants:
    #     txt = "\nSAMPLE: {sample}\n\tNAME: {geneName}\n\tFUNC:{exonicFunc}".format(sample=variant.sample, geneName=variant.geneName, exonicFunc=variant.exonicFunc)
    #     print(txt)


def main(args):
    vcfs = open_file(args.v)
    samples = []
    for vcf in vcfs:
        samples.append(readVCF(vcf))
    showAllSamplesInfo(samples)


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
