#!/bin/bash


function extractLines {
  gawk '\
      BEGIN{OFS="\t"; print "chr\tpos\tref\talt\n----\t----\t----\t----"}
      substr($0, 0, 1)!="#" && \
      /PASS/ \
      {print $1, $2, $4, $5}'
}


function showVariants {
  if [ ${1: -3} == ".gz" ]
    then
      zcat $1 | extractLines
    else
      cat $1 | extractLines
  fi
}

function GC {
  gawk '\
      BEGIN{OFS="\t"; print "\nVar\t#\t%\n----\t----\t----"}
      substr($0, 0, 1)!="#" && \
      /PASS/ \
      {n++; a[$4 ">" $5]++} \
      END{for (pair in a) printf("%s\t%d\t%f\n", pair, a[pair], a[pair]/n)}'
}


function computeCGpourcentage {
  if [ ${1: -3} == ".gz" ]
    then
      zcat $1 | GC
    else
      cat $1 | GC
  fi
}


function countOxog {
  gawk '\
    BEGIN{OFS="\t"} \
    substr($0, 0, 1)!="#" && \
    /PASS/ \
    {n++} \
    {if (($4 == "C" && $5 == "A") || ($4 == "G" && $5 == "T") && $7 == "PASS" && substr($0, 0, 1)!="#") oxog++; print $0} \
    END{printf("\nOxog: %d/%d\n", oxog, n)} \
  '
}


function oxog {
  if [ ${1: -3} == ".gz" ]
    then
      zcat $1 | countOxog
    else
      cat $1 | countOxog
  fi
}

function computeTsTv {
  gawk '\
    BEGIN{OFS="\t"} \
    substr($0, 0, 1)!="#" && \
    /PASS/ \
    {n++} \
    ( ($3=="A" && $4=="G") || ($3=="G" && $4=="A") || ($3=="C" && $4=="T") || ($3=="T" && $4=="C") ) {ts++} END{print "TS/TV:", ts/n  , "\n"}
    '
}



function tstv {
  if [ ${1: -3} == ".gz" ]
    then
      zcat $1 | computeTsTv
    else
      cat $1 | computeTsTv
  fi
}

for var in "$@"
do
  tstv ${var}
  showVariants ${var}
  #computeCGpourcentage ${var}
  #oxog ${var}
done
