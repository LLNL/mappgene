#!/usr/bin/env python3
import os,sys,glob,multiprocessing,time,csv,math
from parsl.app.app import python_app
from os.path import *
from mappgene.subscripts import *

@python_app(executors=['worker'])
def run_ivar(params):
    subject_dir = params['work_dir']
    input_reads = params['input_reads']
    ivar_dir = join(subject_dir, 'ivar')
    raw_dir = join(ivar_dir, 'raw_data')
    smart_remove(raw_dir)
    reads = []

    for input_read in input_reads:
        f = join(raw_dir, basename(input_read))
        smart_copy(input_read, f)
        reads.append(f)

    # Deinterleave if only a single FASTQ was found
    fasta = join(ivar_dir, 'references/PS_1200bp.fasta')
    if len(reads) == 1:
        f = reads[0]
        read1 = replace_extension(f, '_R1.fastq.gz')
        read2 = replace_extension(f, '_R2.fastq.gz')
        stat = replace_extension(f, '.stat')
        run(f"bbduk.sh in={f} out1={read1} out2={read2} ref={fasta} stats={stat} " +
            "k=13 ktrim=l hdist=2 restrictleft=31 statscolumns=5 minlen=65", params)
        # run(f'reformat.sh in={f} out1={read1} out2={read2}', params)
        smart_remove(f)
    elif len(reads) == 2:
        reads.sort()
        read1 = reads[0]
        read2 = reads[1]
    else:
        raise Exception(f'Invalid reads: {reads}')

    output_dir = join(subject_dir, 'ivar_outputs')
    subject = join(output_dir, basename(subject_dir))
    bam = replace_extension(subject, '.bam')
    trimmed = replace_extension(subject, '.trimmed')
    trimmed_sorted = replace_extension(subject, '.trimmed.sorted.bam')
    variants = replace_extension(subject, '.variants')
    masked = replace_extension(subject, '.masked.txt')
    trimmed_masked = replace_extension(subject, '.trimmed.masked.bam')
    final_masked = replace_extension(subject, '.final.masked.variants')
    smart_mkdir(output_dir)
    run(f'bwa index {fasta}', params)
    run(f'bwa mem -t 8 {fasta} {read1} {read2} ' +
        f'| samtools sort -o {bam}', params)
    run(f'ivar trim -b {ivar_dir}/nCoV-2019.scheme.bed -p {trimmed} -i {bam} -e', params)
    run(f'samtools sort {trimmed}.bam -o {trimmed_sorted}', params)

    # call variants with ivar
    run(f'samtools mpileup -aa -A -d 0 -B -Q 0 {trimmed_sorted} | ' +
        f'ivar variants -p {variants} -q 20 -t 0.03 -r {fasta} ' +
        f'-g {ivar_dir}/GCF_009858895.2_ASM985889v3_genomic.gff', params)
    # get primers with mismatches to reference
    run(f'ivar getmasked -i {variants}.tsv -b {ivar_dir}/nCoV-2019.bed ' +
        f'-f {ivar_dir}/nCoV-2019.tsv -p {masked}', params)
    # remove reads with primer mismatches
    run(f'ivar removereads -i {trimmed_sorted} -p {trimmed_masked} ' +
        f'-t {masked} -b {ivar_dir}/nCoV-2019.bed', params)
    # call variants with reads with primer mismatches removed
    run(f'samtools mpileup -aa -A -d 0 -B -Q 0 {trimmed_masked} | ' +
        f'ivar variants -p {final_masked} -q 20 -t 0.03 -r {fasta} ' +
        f'-g {ivar_dir}/GCF_009858895.2_ASM985889v3_genomic.gff', params)





