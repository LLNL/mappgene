#!/usr/bin/env python3
import os,sys,glob,multiprocessing,time,csv,math,pprint
from parsl.app.app import python_app
from os.path import *
from mappgene.subscripts import *

@python_app(executors=['worker'])
def run_vpipe(params):

    subject_dir = params['work_dir']
    subject = basename(subject_dir)
    input_reads = params['input_reads']
    stdout = params['stdout']
    vpipe_dir = join(subject_dir, 'vpipe')
    output_dir = join(subject_dir, 'vpipe_outputs')
    raw_dir = join(vpipe_dir, 'samples/a/b/raw_data')
    smart_remove(raw_dir)
    smart_remove(output_dir)
    smart_mkdir(raw_dir)
    smart_mkdir(output_dir)
    reads = []

    start_time = time.time()
    start_str = f'''
=====================================
Starting V-pipe with subject: {subject}
{get_time_date()}
Arguments: 
{pprint.pformat(params, width=1)}
=====================================
'''
    write(stdout, start_str)
    print(start_str)
    update_permissions(vpipe_dir, params)
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
    # fasta = join(vpipe_dir, 'references/PS_1200bp.fasta')
    fasta = join(vpipe_dir, 'references/NC_045512.2.fasta')
    if len(reads) == 1:
        f = reads[0]
        out1 = replace_extension(f, '_R1.fastq.gz')
        out2 = replace_extension(f, '_R2.fastq.gz')
        stat = replace_extension(f, '.stat')
        run(f"bbduk.sh in={f} out1={out1} out2={out2} ref={fasta} stats={stat} " +
            "k=13 ktrim=l hdist=2 restrictleft=31 statscolumns=5 minlen=65", params)
        smart_remove(f)

    # Update sample.tsv with read length
    smart_remove(join(vpipe_dir, '.snakemake'))
    run(f'cd {vpipe_dir} && ./vpipe --dryrun', params)
    samples_tsv = join(vpipe_dir, 'samples.tsv')
    with open(samples_tsv, 'r') as fr, open(f'{samples_tsv}.tmp', 'w') as fw:
        r = csv.reader(fr, delimiter='\t')
        w = csv.writer(fw, delimiter='\t')
        for row in r:
            row.append(params['read_length'])
            w.writerow(row)
    smart_copy(f'{samples_tsv}.tmp', samples_tsv)

    # Run V-pipe analysis
    ncores = int(math.floor(multiprocessing.cpu_count() / 2))
    run(f'cd {vpipe_dir} && ./vpipe --cores {ncores} --use-conda --unlock', params)
    run(f'cd {vpipe_dir} && ./vpipe --cores {ncores} --use-conda', params)
    time.sleep(10)

    # Run snpEff postprocessing
    vcf_s0 = join(vpipe_dir, 'samples/a/b/variants/SNVs/snvs.vcf')
    vcf_s1 = join(output_dir, f'{subject}.vpipe.vcf')
    vcf_s2 = join(output_dir, f'{subject}.vpipe.snpEFF.vcf')
    vcf_s3 = join(output_dir, f'{subject}.vpipe.snpSIFT.txt')
    run(f'sed "s/MN908947.3/NC_045512.2/g" {vcf_s0} > {vcf_s1}', params)
    run(f'java -Xmx8g -jar /opt/snpEff/snpEff.jar NC_045512.2 {vcf_s1} > {vcf_s2}', params)
    run(f'cat {vcf_s2} | /opt/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /opt/snpEff/SnpSift.jar ' +
        f' extractFields - CHROM POS REF ALT AF DP "ANN[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].EFFECT" ' +
        f' "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].AA_POS" "ANN[*].GENE" > {vcf_s3}', params)

    # Copy outputs to accessible location
    smart_copy(join(vpipe_dir, 'samples/a/b/alignments'), join(output_dir, 'alignments'))
    smart_copy(join(vpipe_dir, 'samples/a/b/variants'), join(output_dir, 'variants'))
    smart_copy(join(vpipe_dir, 'samples/a/b/extracted_data'), join(output_dir, 'extracted_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(vpipe_dir, 'samples/a/b/preprocessed_data'), join(output_dir, 'preprocessed_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(vpipe_dir, 'samples/a/b/raw_data'), join(output_dir, 'raw_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(vpipe_dir, 'samples/a/b/references'), join(output_dir, 'references'), exclude=['*.gz', '*.fasta'])

    # Clear extra files
    smart_remove('snpEff_genes.txt')
    smart_remove('snpEff_summary.html')

    update_permissions(vpipe_dir, params)
    update_permissions(output_dir, params)
    finish_str = f'''
=====================================
Finished V-pipe with subject: {subject}
{get_time_date()}
Arguments: 
{pprint.pformat(params, width=1)}
Total time: {get_time_string(time.time() - start_time)} (HH:MM:SS)
=====================================
'''
    write(stdout, finish_str)
    print(finish_str)

