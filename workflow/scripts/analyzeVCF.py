import re
import argparse


class Filter():
    all_chrom = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "Y"]
    def __init__(self, samples=[], chrom=all_chrom, germAF = 0, diagAF = 0.05, tlod = 8, roq = 0, f1r2 = 3, f2r1 = 3, mbq = 20, mmq = 20):
        self.samples = samples
        self.chrom = chrom
        self.germAF = germAF
        self.diagAF = diagAF
        self.tlod  = tlod
        self.roq = roq
        self.f1r2 = f1r2
        self.f2r1 = f2r1
        self.mbq = mbq
        self.mmq = mmq


    def applyFilter(self, samples):
        print('#SAMPLE\tCHROM\tPOS\tREF\tALT\tGENE_NAME\tFUNC\tGERM_AF\tDIAG_AF\tTLOD\tNLOD\tAD_GERM\tROQ\tF1R2\tF2R1\tMBQ\tMMQ')
        for sample in samples:
            if sample in self.samples:
                continue
            else:
                for variant in samples[sample]:
                    variant_pass = True
                    if not variant.chr in self.chrom:
                        variant_pass = False
                    if float(variant.somaticAF) < float(self.diagAF):
                        variant_pass = False
                        #print("Somatic VAF")
                    if float(variant.germlineAF) < float(self.germAF):
                        variant_pass = False
                        #print("Germline VAF")
                    if float(variant.TLOD) < float(self.tlod):
                        variant_pass = False
                        #print("TLOD")
                    if float(variant.ROQ) < float(self.roq):
                        variant_pass = False
                        #print("ROQ")
                    if int(variant.F1R2) < int(self.f1r2):
                        variant_pass = False
                        #print("F1R2")
                    if int(variant.F2R1) < int(self.f2r1):
                        variant_pass = False
                        #print("F2R1")
                    if int(variant.MMQ) < int(self.mmq):
                        variant_pass = False
                        #print("MMQ")
                    if int(variant.MBQ) < int(self.mbq):
                        variant_pass = False
                        #print("MBQ")
                    if variant_pass == True:
                        txt = "{sample_name}\t{chr}\t{pos}\t{ref}\t{alt}\t{geneName}\t{exonicFunc}\t{germAF}\t{diagAF}\t{tLOD}\t{nLOD}\t{ad}\t{roq}\t{f1r2}\t{f2r1}\t{mbq}\t{mmq}".format(mmq=variant.MMQ, mbq=variant.MBQ, f2r1=variant.F2R1, f1r2=variant.F1R2, roq=variant.ROQ, ad=variant.AD, nLOD=variant.NLOD, tLOD=variant.TLOD, germAF=variant.germlineAF, diagAF=variant.somaticAF, sample_name=sample, geneName=variant.geneName, exonicFunc=variant.exonicFunc, chr=variant.chr, pos=variant.pos, ref=variant.ref, alt=variant.alt)
                        print(txt)


class Sample():
    def __init__(self, sample_name, variants):

        CLONAL_INTERFERENCE_GENES = ["KIT", "FLT3", "NRAS", "KRAS", "JAK2", "CBL"]  # temporaire

        NOT_CAUSALS_VARIATIONS = ["frameshift_deletion", "frameshift_insertion", "None", "synonymous_SNV", "unknown"]

        def getClonalInterferenceVariants(variants):
            cpt  = 0
            genes = []
            nTotal = 0
            for variant in variants:
                nTotal += 1
                if (variant.geneName in CLONAL_INTERFERENCE_GENES) and (variant.exonicFunc not in NOT_CAUSALS_VARIATIONS):
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


        def getAlleleFrenquencies(line):
            d, g = line.split("\t")[9:11]
            return [d.split(":")[2], g.split(":")[2]]


        def getAD(line):
            d = line.split("\t")[9]
            return d.split(":")[1].split(",")[1]


        def getTLOD(line):
            m = re.search("TLOD=([0-9]+(\.[0-9][0-9]?)?)", line)
            if m:
                return m.group(1)

        def getNLOD(line):
            m = re.search("NLOD=([0-9]+(\.[0-9][0-9]?)?);", line)
            if m:
                return m.group(1)


        def getROQ(line):
            m = re.search("ROQ=([0-9]+(\.[0-9][0-9]?)?);", line)
            if m:
                return m.group(1)


        def getF1R2(line):
            d = line.split("\t")[9]
            return d.split(":")[4].split(",")[1]


        def getF2R1(line):
            d = line.split("\t")[9]
            return d.split(":")[5].split(",")[1]


        def getMBQ(line):
            m = re.search("MBQ=([0-9]+(\.[0-9][0-9]?)?),([0-9]+(\.[0-9][0-9]?)?);", line)
            if m:
                return m.group(3)


        def getMMQ(line):
            m = re.search("MMQ=([0-9]+(\.[0-9][0-9]?)?),([0-9]+(\.[0-9][0-9]?)?);", line)
            if m:
                return m.group(3)



        self.sample = get_sample(vcf_file)
        self.chr = get_chr(line)
        self.pos = get_pos(line)
        self.ref = get_ref(line)
        self.alt = get_alt(line)
        self.geneName = get_geneName(line)
        self.exonicFunc = get_exonicFunc(line)
        self.somaticAF, self.germlineAF = getAlleleFrenquencies(line)
        self.TLOD = getTLOD(line)
        self.NLOD = getNLOD(line)
        self.AD = getAD(line)
        self.ROQ = getROQ(line)
        self.F1R2 = getF1R2(line)
        self.F2R1 = getF2R1(line)
        self.MBQ = getMBQ(line)
        self.MMQ = getMMQ(line)


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
    print('#SAMPLE\tCHROM\tPOS\tREF\tALT\tGENE_NAME\tFUNC\tGERM_AF\tDIAG_AF\tTLOD\tNLOD\tAD_GERM\tROQ\tF1R2\tF2R1\tMBQ\tMMQ\n')
    for sample in samples:
        for variant in samples[sample]:
            txt = "{sample_name}\t{chr}\t{pos}\t{ref}\t{alt}\t{geneName}\t{exonicFunc}\t{germAF}\t{diagAF}\t{tLOD}\t{nLOD}\t{ad}\t{roq}\t{f1r2}\t{f2r1}\t{mbq}\t{mmq}\n".format(mmq=variant.MMQ, mbq=variant.MBQ, f2r1=variant.F2R1, f1r2=variant.F1R2, roq=variant.ROQ, ad=variant.AD, nLOD=variant.NLOD, tLOD=variant.TLOD, germAF=variant.germlineAF, diagAF=variant.somaticAF, sample_name=sample, geneName=variant.geneName, exonicFunc=variant.exonicFunc, chr=variant.chr, pos=variant.pos, ref=variant.ref, alt=variant.alt)
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
    if args.g:
        grouped_samples = assignCategorie(samples)
        showGroupedSamples(grouped_samples)
    if args.f:
        filter = Filter()
        filter.applyFilter(samples)
    else:
        showAllSamplesInfo(samples)


if __name__== '__main__':
    parser = argparse.ArgumentParser(description='Analyze VCF.')
    parser.add_argument('-v', default=None, required=True, type=str, help="File with all VCFs name.")
    parser.add_argument('-g', default=False, required=False, type=bool, help="Determine and show sample's group.")
    # parser.add_argument('-e', default='SJCBF', type=str, help="Experience. Default: SJCBF.")
    # parser.add_argument('-j', default=True, type=bool, help="Create JSON. Default: True")
    # parser.add_argument('-p', default='/data1/scratch/pamesl/projet_cbf/data/bam/', type=str, help="Path to bam files repository.")
    # parser.add_argument('-l', default=1, type=int, help="Number of files to download. Default: 1.")
    # parser.add_argument('-t', default=None, type=str, help="Type of files to download. For exemple: 'D' or 'G'.", required=True)
    parser.add_argument('-f', default=False, type=bool, help="Filter vcf.", required=False)

    args = parser.parse_args()

    main(args)
