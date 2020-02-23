import argparse

class Variant():
    def __init__(self, line, vcf_file):

        def get_chr(line):
            return line.split("\t")[0]


        def get_pos(line):
            return line.split("\t")[1]


        def get_ref(line):
            return line.split("\t")[3]


        def get_alt(line):
            return line.split("\t")[4]


        def get_variantCaller(vcf_file):
            return vcf_file.split("/")[3].split("_")[0]

        self.chr = get_chr(line)
        self.pos = get_pos(line)
        self.ref = get_ref(line)
        self.alt = get_alt(line)
        self.vc = get_variantCaller(vcf_file)


def open_file(file_path):
    if file_path[-1] == "\n":
        file = open(file_path[0:-1], 'r')
    else:
        file = open(file_path, 'r')
    return file.readlines()


def readVCF(vcf_file):
    lines = open_file(vcf_file)
    return [Variant(line, vcf_file) for line in lines if line[0] != "#"]


def printHeader():
    print("CHR\tPOS\tREF\tALT\tVC")

def printVariant(v):
    txt = "{chr}\t{pos}\t{ref}\t{alt}\t{vc}"
    print(txt.format(chr=v.chr, pos=v.pos, ref=v.ref, alt=v.alt, vc=v.vc))

def main(args):
    vcfs = open_file(args.v)
    variants = []
    for vcf in vcfs:
        variants.extend(readVCF(vcf))
    #printHeader()
    for variant in variants:
        printVariant(variant)



if __name__== '__main__':
    parser = argparse.ArgumentParser(description='Analyze VCF.')
    parser.add_argument('-v', default=None, required=True, type=str, help="File with all VCFs name.")

    args = parser.parse_args()

    main(args)
