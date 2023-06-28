#!/usr/bin/env bash

if ! which makeblastdb > /dev/null; then
    echo "Needs NCBI BLAST+ >= 2.2.30"
    exit 1
fi

# if ! parallel --version | grep -q -s GNU ; then
#     echo "Needs GNU parallel"
#     exit 2
# fi

# parallel --skip-first-line --colsep $'\t' \
#     "echo -e '>{1}\n{2}\n'" \
#     < Primers\ for\ amplicon-based\ sequencing\ -\ ZIKV\ -\ Asia_America\ -\ 400bp.tsv \
#     > zibra.fa

map_primers(){
    local _ref="$1"
    local _primers="$2"

    makeblastdb -in "$_ref" -input_type fasta -dbtype nucl > /dev/null

    # blast the primers, using strict parameters
    #
    # output is pretty much in BED format, it needs a column for the pool number(?)
    # to match the other files we have for ivar input. however, i don't think the
    # pool number nor the seq is used by ivar.
    blastn -query "$_primers" -db "$_ref" \
        -word_size 8 -qcov_hsp_perc 0.8 -max_hsps 1 -ungapped \
        -outfmt "6 sseqid sstart send qseqid sstrand qseq qlen length evalue pident qcovs"
}

map_zibra_primers(){
    local _acc="$1"
    local _primers_dir="primers_zibra_${_acc}bp"

    mkdir -p "$_primers_dir"
    map_primers "references/${_acc}.fasta" zibra.fa > "${_primers_dir}/blast.tsv"

    # drop score columns and add dummy pool numbers (0)
    cd "$_primers_dir"
    awk -v OFS=$'\t' '{
                          gsub(/plus/, "+", $5);
                          gsub(/minus/, "-", $5);
                          print $1, $2, $3, $4, "0", $5, $6
                      }' < blast.tsv \
    | sort -nk2 > primers.bed
    sort -Vk4 primers.bed \
    | awk '  NR%2  {prev=$4}
           !(NR%2) {print prev, "\t", $4}' > pairinfo.tsv
    cd -
}

map_zibra_primers KU501215.1
map_zibra_primers KU955593.1
map_zibra_primers KX087101.3
