#!/usr/bin/env python3
from parsl.app.app import python_app

@python_app(executors=['worker'], cache=True)
def run_vpipe(params, inputs=[]):
    import math,multiprocessing,time,csv,os,sys
    from os.path import basename,join

    sys.path.append(subscripts_dir)
    from subscripts.utilities import smart_copy,smart_mkdir,smart_remove,run

    # Setup parameters and input data
    work_dir = join(output_dir, 'work_dir')
    smart_mkdir(work_dir)
    params['sdir'] = work_dir
    params['stdout'] = join(work_dir, 'worker.stdout')
    smart_copy(params['git_dir'], work_dir)
    work_input_dir = join(work_dir, 'samples/a/b/raw_data')

    # Run fixq.sh
    for f in inputs:
        work_f = join(work_input_dir, basename(f))
        work_f2 = join(work_input_dir, 'tmp_' + basename(f))
        smart_copy(f, work_f)
        run(f'zcat {work_f} | awk \'NR%4 == 0 {{ gsub(\\"F\\", \\"?\\"); gsub(\\":\\", \\"5\\") }}1\'' +
            f' | gzip -c > {work_f2}', params)
        smart_remove(work_f)
        os.rename(work_f2, work_f)

    # Deinterleave if only a single FASTQ was found
    if len(inputs) == 1:
        f = inputs[0]
        work_f = join(work_input_dir, basename(f))
        work_r1 = join(work_input_dir, basename(f).replace('.fastq.gz', '_R1.fastq.gz'))
        work_r2 = join(work_input_dir, basename(f).replace('.fastq.gz', '_R2.fastq.gz'))
        work_stat = join(work_input_dir, basename(f).replace('.fastq.gz', '.stat'))
        fasta = join(work_dir, 'references/PS_1200bp.fasta')
        run(f"bbduk.sh in={work_f} out1={work_r1} out2={work_r2} ref={fasta} stats={work_stat} " +
            "k=13 ktrim=l hdist=2 restrictleft=31 statscolumns=5 minlen=65", params)
        smart_remove(work_f)

    # Update sample.tsv with read length
    run(f'cd {work_dir} && ./vpipe --dryrun', params)
    samples_tsv = join(work_dir, 'samples.tsv')
    with open(samples_tsv, 'r') as fr, open(f'{samples_tsv}.tmp', 'w') as fw:
        r = csv.reader(fr, delimiter='\t')
        w = csv.writer(fw, delimiter='\t')
        for row in r:
            row.append(params['read_length'])
            w.writerow(row)
    smart_copy(f'{samples_tsv}.tmp', samples_tsv)

    # Run V-pipe analysis
    ncores = int(math.floor(multiprocessing.cpu_count() / 2))
    run(f'cd {work_dir} && ./vpipe --cores {ncores} --use-conda', params)
    time.sleep(10)

    smart_copy('/usr/WS2/moon15/old/mappgene/output/example/work_dir/samples',
        '/usr/WS2/moon15/mappgene/output/example/work_dir/samples')
    vcf_s0 = join(work_dir, 'samples/a/b/variants/SNVs/snvs.vcf')
    vcf_s1 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.vcf')
    vcf_s2 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.2.snpEFF.vcf')
    vcf_s3 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.2.snpSIFT.txt')

    run(f'sed "s/MN908947.3/NC_045512.2/g" {vcf_s0} > {vcf_s1}', params)
    # run('java -jar /opt/snpEff/snpEff.jar download -v NC_045512.2', params)
    run(f'java -Xmx8g -jar /opt/snpEff/snpEff.jar NC_045512.2 {vcf_s1} > {vcf_s2}', params)
    run(f'cat {vcf_s2} | /opt/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /opt/snpEff/SnpSift.jar ' +
        f' extractFields - CHROM POS REF ALT AF DP "ANN[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].EFFECT" ' +
        f' "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].AA_POS" "ANN[*].GENE" > {vcf_s3}', params)

    smart_copy(join(work_dir, 'samples/a/b/alignments'), join(output_dir, 'alignments'))
    smart_copy(join(work_dir, 'samples/a/b/variants'), join(output_dir, 'variants'))
    smart_copy(join(work_dir, 'samples/a/b/extracted_data'), join(output_dir, 'extracted_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(work_dir, 'samples/a/b/preprocessed_data'), join(output_dir, 'preprocessed_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(work_dir, 'samples/a/b/raw_data'), join(output_dir, 'raw_data'), exclude=['*.gz', '*.fasta'])
    smart_copy(join(work_dir, 'samples/a/b/references'), join(output_dir, 'references'), exclude=['*.gz', '*.fasta'])

