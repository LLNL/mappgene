#!/usr/bin/env python3
import os,sys,glob,multiprocessing,time,csv,math,pprint
from parsl.app.app import python_app
from os.path import *
from mappgene.subscripts import *

@python_app(executors=['worker'])
def run_ivar(params):

    subject_dir = params['work_dir']
    subject = basename(subject_dir)
    input_reads = params['input_reads']
    variant_frequency = params['variant_frequency']
    stdout = params['stdout']
    ivar_dir = join(subject_dir, 'ivar')
    output_dir = join(subject_dir, 'ivar_outputs')
    alignments_dir = join(output_dir, 'alignments')
    raw_dir = join(ivar_dir, 'raw_data')
    smart_remove(raw_dir)
    smart_remove(output_dir)
    smart_mkdir(raw_dir)
    smart_mkdir(output_dir)
    smart_mkdir(alignments_dir)
    reads = []

    start_time = time.time()
    start_str = f'''
=====================================
Starting iVar with subject: {subject}
{get_time_date()}
Arguments: 
{pprint.pformat(params, width=1)}
=====================================
'''
    write(stdout, start_str)
    print(start_str)
    update_permissions(ivar_dir, params)
    update_permissions(output_dir, params)

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
        run(f'reformat.sh in={f} out1={read1} out2={read2}', params)
        smart_remove(f)
    elif len(reads) == 2:
        reads.sort()
        read1 = reads[0]
        read2 = reads[1]
    else:
        raise Exception(f'Invalid reads: {reads}')

    
    align_prefix = join(alignments_dir, subject)
    bam = replace_extension(align_prefix, '.bam')
    trimmed = replace_extension(align_prefix, '.trimmed')
    trimmed_sorted = replace_extension(align_prefix, '.trimmed.sorted.bam')
    variants = replace_extension(align_prefix, '.variants')
    noinsertions = replace_extension(align_prefix, '.noins.variants')
    masked = replace_extension(align_prefix, '.masked.txt')
    trimmed_masked = replace_extension(align_prefix, '.trimmed.masked.bam')
    final_masked = replace_extension(align_prefix, '.final.masked.variants')
    lofreq_bam = replace_extension(align_prefix, '.lofreq.bam')
    vcf_s0 = replace_extension(align_prefix, '.vcf')
    tsv = replace_extension(align_prefix, '.final.masked.variants.tsv')
    output_tsv = join(output_dir, f'{subject}.ivar.final.masked.variants.tsv')
    run(f'bwa index {fasta}', params)
    run(f'bwa mem -t 8 {fasta} {read1} {read2} | samtools sort -o {bam}', params)
    run(f'ivar trim -b {ivar_dir}/nCoV-2019.scheme.bed -p {trimmed} -i {bam} -e', params)
    run(f'samtools sort {trimmed}.bam -o {trimmed_sorted}', params)

    # call variants with ivar (produces {subject}.variants.tsv)
    run(f'samtools mpileup -aa -A -d 0 -B -Q 0 {trimmed_sorted} | ' +
        f'ivar variants -p {variants} -q 20 -t {variant_frequency} -r {fasta} ' +
        f'-g {ivar_dir}/GCF_009858895.2_ASM985889v3_genomic.gff', params)
    
    # remove low quality insertions because we want to ignore most mismatches
    # to primers that are insertions (produces {subject}.noins.variants.tsv)
    run(f"awk \'! (\\$4 ~ /^\\+/ && \\$10 >= 20) {{ print }}\' < {variants}.tsv > {noinsertions}.tsv", params)
    
    # get primers with mismatches to reference (produces {subject}.masked.txt)
    run(f'ivar getmasked -i {noinsertions}.tsv -b {ivar_dir}/nCoV-2019.bed ' +
        f'-f {ivar_dir}/nCoV-2019.tsv -p {masked}', params)
    
    # remove reads with primer mismatches (produces {subject}.trimmed.masked.bam)
    run(f'ivar removereads -i {trimmed_sorted} -p {trimmed_masked} ' +
        f'-t {masked} -b {ivar_dir}/nCoV-2019.bed', params)
    
    # call variants with reads with primer mismatches removed (produces {subject}.final.masked.variants)
    run(f'samtools mpileup -aa -A -d 0 -B -Q 0 {trimmed_masked} | ' +
        f'ivar variants -p {final_masked} -q 20 -t {variant_frequency} -r {fasta} ' +
        f'-g {ivar_dir}/GCF_009858895.2_ASM985889v3_genomic.gff', params)
    smart_copy(tsv, output_tsv)

    # # use lofreq to convert bam to vcf
    run(f'lofreq indelqual --dindel -f {fasta} -o {lofreq_bam} --verbose {trimmed_masked}', params)
    run(f'samtools index {lofreq_bam}', params)
    run(f'lofreq call --call-indels -f {fasta} -o {vcf_s0} --verbose {lofreq_bam}', params)

    # Run snpEff postprocessing
    vcf_s1 = join(output_dir, f'{subject}.ivar.vcf')
    vcf_s2 = join(output_dir, f'{subject}.ivar.snpEFF.vcf')
    vcf_s3 = join(output_dir, f'{subject}.ivar.snpSIFT.txt')
    run(f'sed "s/MN908947.3/NC_045512.2/g" {vcf_s0} > {vcf_s1}', params)
    run(f'java -Xmx8g -jar /opt/snpEff/snpEff.jar NC_045512.2 {vcf_s1} > {vcf_s2}', params)
    run(f'cat {vcf_s2} | /opt/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /opt/snpEff/SnpSift.jar ' +
        f' extractFields - CHROM POS REF ALT AF DP "ANN[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].EFFECT" ' +
        f' "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].AA_POS" "ANN[*].GENE" > {vcf_s3}', params)

    # Clear extra files
    smart_remove('snpEff_genes.txt')
    smart_remove('snpEff_summary.html')

    update_permissions(ivar_dir, params)
    update_permissions(output_dir, params)
    finish_str = f'''
=====================================
Finished iVar with subject: {subject}
{get_time_date()}
Arguments: 
{pprint.pformat(params, width=1)}
Total time: {get_time_string(time.time() - start_time)} (HH:MM:SS)
=====================================
'''
    write(stdout, finish_str)
    print(finish_str)
