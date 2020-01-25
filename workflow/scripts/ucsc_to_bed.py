#  UCSC To BED

import csv
import sys
import re


with open(sys.argv[1]) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        chr = row[1]
        starts = row[8]
        ends = row[9]
        if re.match('^chr[XY12345678910M]+$', chr) is not None:
            ends = ends.split(',')
            starts = starts.split(',')
            length = len(ends)
            if length == len(starts):
                for i in range(0, length-1):
                    print("{chr}\t{start}\t{end}".format(chr=chr, start=starts[i], end=ends[i]))
