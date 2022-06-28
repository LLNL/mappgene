#!/bin/sh

cp -v SARS-CoV-2.scheme.bed nCoV-2019.scheme.bed
cp -v SARS-CoV-2.primer.bed nCoV-2019.bed

echo creating nCoV-2019.tsv from SARS-CoV-2.scheme.bed
echo assuming odd lines are LEFT, even lines are RIGHT
awk 'NR%2 {prev=$4} \
     !(NR%2) {print prev, "\t", $4}' < SARS-CoV-2.scheme.bed > nCoV-2019.tsv
