#!/bin/bash

awk '($10 == "2") {print $1"\t"$2"\t"$3}' $1 | tail -n +2 > $2
