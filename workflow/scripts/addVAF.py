
import sys

index = {"A":4, "C":5, "G":6, "T":7}

for line in sys.stdin:
    splited = line.split("\t")
    ref = splited[3]
    alt = splited[4]
    refCount = splited[9].split(":")[index[ref]]
    altCount = splited[10].split(":")[index[alt]]
    tier1RefCounts = refCount.split(",")[0]
    tier1AltCounts = altCount.split(",")[0]
    vaf = float(tier1AltCounts) / ((float(tier1AltCounts) + float(tier1RefCounts)))
    print(splited[0] + "\t" + splited[1] + "\t" + splited[3] + "\t" + splited[4] + "\t" + str(vaf))
