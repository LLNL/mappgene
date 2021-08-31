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
    output_dir = join(subject_dir, 'ivar_outputs')
    alignments_dir = join(output_dir, 'alignments')
    raw_dir = join(ivar_dir, 'raw_data')
    smart_remove(raw_dir)
    smart_remove(output_dir)
    smart_mkdir(output_dir)
    smart_mkdir(alignments_dir)
    reads = []

    # Run fixq.sh
    for input_read in input_reads:
        tmp_f = join(raw_dir, 'tmp_' + basename(input_read))
        f = join(raw_dir, basename(input_read))
        smart_copy(input_read, f)
        run(f'zcat {f} | awk \'NR%4 == 0 {{ gsub(\\"F\\", \\"?\\"); gsub(\\":\\", \\"5\\") }}1\'' +
            f' | gzip -c > {tmp_f}', params)
        if exists(tmp_f):
            smart_remove(f)
            os.rename(tmp_f, f)
        reads.append(f)

    # Deinterleave if only a single FASTQ was found
    # fasta = join(ivar_dir, 'references/PS_1200bp.fasta')
    fasta = join(ivar_dir, 'references/NC_045512.2.fasta')
    if len(reads) == 1:
        f = reads[0]
        read1 = replace_extension(f, '_R1.fastq.gz')
        read2 = replace_extension(f, '_R2.fastq.gz')
        stat = replace_extension(f, '.stat')
        run(f"bbduk.sh in={f} out1={read1} out2={read2} ref={fasta} stats={stat} " +
            "k=13 ktrim=l hdist=2 restrictleft=31 statscolumns=5 minlen=65", params)
        smart_remove(f)
    elif len(reads) == 2:
        reads.sort()
        read1 = reads[0]
        read2 = reads[1]
    else:
        raise Exception(f'Invalid reads: {reads}')

    
    subject = join(alignments_dir, basename(subject_dir))
    bam = replace_extension(subject, '.bam')
    trimmed = replace_extension(subject, '.trimmed')
    trimmed_sorted = replace_extension(subject, '.trimmed.sorted.bam')
    variants = replace_extension(subject, '.variants')
    masked = replace_extension(subject, '.masked.txt')
    trimmed_masked = replace_extension(subject, '.trimmed.masked.bam')
    final_masked = replace_extension(subject, '.final.masked.variants')
    lofreq_bam = replace_extension(subject, '.lofreq.bam')
    vcf_s0 = replace_extension(subject, '.vcf')
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

    # # use lofreq to convert bam to vcf
    run(f'lofreq indelqual --dindel -f {fasta} -o {lofreq_bam} --verbose {trimmed_masked}', params)
    run(f'samtools index {lofreq_bam}', params)
    run(f'lofreq call --call-indels -f {fasta} -o {vcf_s0} --verbose {lofreq_bam}', params)

    # Run snpEff postprocessing
    vcf_s1 = join(output_dir, 'snvs_NC_045512.vcf')
    vcf_s2 = join(output_dir, 'snvs_NC_045512.2.snpEFF.vcf')
    vcf_s3 = join(output_dir, 'snvs_NC_045512.2.snpSIFT.txt')
    run(f'sed "s/MN908947.3/NC_045512.2/g" {vcf_s0} > {vcf_s1}', params)
    run(f'java -Xmx8g -jar /opt/snpEff/snpEff.jar NC_045512.2 {vcf_s1} > {vcf_s2}', params)
    run(f'cat {vcf_s2} | /opt/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /opt/snpEff/SnpSift.jar ' +
        f' extractFields - CHROM POS REF ALT AF DP "ANN[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].EFFECT" ' +
        f' "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].AA_POS" "ANN[*].GENE" > {vcf_s3}', params)

    # Clear extra files
    smart_remove('snpEff_genes.txt')
    smart_remove('snpEff_summary.html')

