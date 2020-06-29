#!/bin/bash


function cutVCF {
  if [ ${1: -3} == ".gz" ]
    then
      zcat $1 | gawk '\
      {OFS="\t"}
      substr($0, 0, 1)!="#" && \
      /PASS/ \
      {print $1, $2, $4, $5}'
    else
      cat $1 | gawk '\
      {OFS="\t"}
      substr($0, 0, 1)!="#" && \
      /PASS/ \
      {print $1, $2, $4, $5}'
  fi
}


for var in "$@"
do
  cutVCF ${var}
done | sort -n | uniq -c | awk '{OFS="\t"} ($1==2) {print $2, $3, $4, $5, $6}'
