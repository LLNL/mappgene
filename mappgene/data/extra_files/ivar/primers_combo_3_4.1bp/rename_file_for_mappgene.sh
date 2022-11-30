#!/bin/sh
# concatenates v3 and v4 primers

_V3_DIR=../primers_400bp
_V4_1_DIR=../primers_v4.1bp

cat $_V3_DIR/nCoV-2019.scheme.bed $_V4_1_DIR/nCoV-2019.scheme.bed \
    | sort -k 1,1 -k2,2n > nCoV-2019.scheme.bed
cat $_V3_DIR/nCoV-2019.bed        $_V4_1_DIR/nCoV-2019.bed \
    | sort -k 1,1 -k2,2n > nCoV-2019.bed
awk '  NR%2  {prev=$4} !(NR%2) {print prev, "\t", $4}' $_V3_DIR/nCoV-2019.scheme.bed \
    | cat - $_V4_1_DIR/nCoV-2019.tsv > nCoV-2019.tsv

