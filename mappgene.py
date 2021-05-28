#!/usr/bin/env python3
import argparse,multiprocessing,parsl,getpass,socket,json,sys,re,glob,math
from distutils.dir_util import copy_tree
from parsl.app.app import python_app,bash_app
from parsl.executors import ThreadPoolExecutor,HighThroughputExecutor
from parsl.providers import LocalProvider,SlurmProvider
from parsl.addresses import address_by_hostname,address_by_route
from os.path import exists,join,split,splitext,abspath,basename,islink,isdir
from subscripts.utilities import *

class ArgsObject:
    def __init__(self, **entries):
        self.__dict__.update(entries)

if len(sys.argv) == 2 and sys.argv[1] not in ['-h', '--help', '-f', '--force', '--local']:
    config_json = sys.argv[1]
    with open(config_json) as f:
        raw_args = json.load(f)
    args = ArgsObject(**raw_args)
else:
    parser = argparse.ArgumentParser(description='Generate connectome and edge density images',
        usage=""" """)

    parser.add_argument('--input_dirs', '-i', help='Path to inputs. Subdirectories must contain subject genomes in FASTQ format.')
    parser.add_argument('--output_dirs', '-o', help='Path to outputs.')
    parser.add_argument('--nnodes', '-n', help='Number of nodes.')
    parser.add_argument('--local', help='Run only on local threads (for debugging).', action='store_true')
    parser.add_argument('--bank', '-b', help='Bank to charge for jobs.')
    parser.add_argument('--partition', '-p', help='Scheduler partition to assign jobs.')
    parser.add_argument('--force', '-f', help='Overwrite existing outputs.', action='store_true')
    parser.add_argument('--container', help='Path to Singularity container image.')
    parser.add_argument('--walltime', '-t', help='Walltime in format HH:MM:SS.')
    parser.add_argument('--single_subject', help='Run with just a single subject from the input directory.')
    parser.add_argument('--read_length', help='Read length in sample.tsv (see cbg-ethz.github.io/V-pipe/tutorial/sars-cov2).')
    args = parser.parse_args()

pending_args = args.__dict__.copy()
parse_default('input_dirs', 'input/', args, pending_args)
parse_default('output_dirs', 'output/', args, pending_args)
parse_default('bank', 'asccasc', args, pending_args)
parse_default('partition', 'pbatch', args, pending_args)
parse_default('force', False, args, pending_args)
parse_default('local', False, args, pending_args)
parse_default('container', "container/image.sif", args, pending_args)
parse_default('nnodes', 1, args, pending_args)
parse_default('walltime', '11:59:00', args, pending_args)
parse_default('single_subject', '', args, pending_args)
parse_default('read_length', 250, args, pending_args)

if __name__ == '__main__':

    # Setup V-pipe repo
    tmp_dir = abspath('tmp')
    smart_remove(tmp_dir)
    smart_mkdir(tmp_dir)
    git_dir = join(tmp_dir, 'vpipe')
    git_params = {'sdir':git_dir, 'container':args.container}
    # run(f'git clone https://github.com/cbg-ethz/V-pipe.git {git_dir}', git_params)
    run(f'cp -rf /opt/vpipe {git_dir}', git_params)
    copy_tree('vpipe_files', join(git_dir))
    run(f'sh -c "cd {git_dir} && sh init_project.sh" || true', git_params)

    if args.local:
        executor = ThreadPoolExecutor(label="worker")
    else:
        executor = HighThroughputExecutor(
            label="worker",
            address=address_by_hostname(),
            provider=SlurmProvider(
                args.partition,
                launcher=parsl.launchers.SrunLauncher(),
                nodes_per_block=int(args.nnodes),
                init_blocks=1,
                max_blocks=1,
                worker_init=f"export PYTHONPATH=$PYTHONPATH:{os.getcwd()}",
                walltime=args.walltime,
                scheduler_options="#SBATCH --exclusive\n#SBATCH -A {}\n".format(args.bank),
                move_files=False,
            ),
        )
    params = {
        'container': abspath(args.container),
        'git_dir': abspath(git_dir),
        'read_length': args.read_length,
    }

    config = parsl.config.Config(executors=[executor])
    parsl.set_stream_logger()
    parsl.load(config)


    @python_app(executors=['worker'], cache=True)
    def run_worker(input_dir, output_dir, params):
        import math,multiprocessing,glob,time,csv
        from os.path import basename,join
        from subscripts.utilities import smart_copy,smart_mkdir,run

        # Setup parameters and input data
        work_dir = join(output_dir, 'work_dir')
        smart_mkdir(work_dir)
        params['sdir'] = work_dir
        params['stdout'] = join(work_dir, 'worker.stdout')
        smart_copy(params['git_dir'], work_dir)
        inputs = glob.glob(join(input_dir, '*.fastq.gz'))
        work_input_dir = join(work_dir, 'samples/a/b/raw_data')

        # # Deinterleave if only a single FASTQ is found
        if len(inputs) == 1:
            f = inputs[0]
            work_f = join(work_input_dir, basename(f))
            work_r1 = join(work_input_dir, basename(f).replace('.fastq.gz', '_R1.fastq.gz'))
            work_r2 = join(work_input_dir, basename(f).replace('.fastq.gz', '_R2.fastq.gz'))
            work_stat = join(work_input_dir, basename(f).replace('.fastq.gz', '.stat'))
            fasta = join(work_dir, 'references/PS_1200bp.fasta')
            smart_copy(f, work_f)
            run(f"bbduk.sh in={work_f} out1={work_r1} out2={work_r2} ref={fasta} stats={work_stat} " +
                "k=13 ktrim=l hdist=0 restrictleft=31 statscolumns=5", params)
            smart_remove(work_f)
        else:
            for f in inputs:
                smart_copy(f, join(work_input_dir, basename(f)))
        return

        # Update sample.tsv with read length
        run(f'sh -c "cd {work_dir} && ./vpipe --dryrun"', params)
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
        run(f'sh -c "cd {work_dir} && ./vpipe --cores {ncores} --use-conda"', params)
        time.sleep(10)

        smart_copy('/usr/WS2/moon15/old/mappgene/output/example/work_dir/samples',
            '/usr/WS2/moon15/mappgene/output/example/work_dir/samples')
        vcf_s0 = join(work_dir, 'samples/a/b/variants/SNVs/snvs.vcf')
        vcf_s1 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.vcf')
        vcf_s2 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.2.snpEFF.vcf')
        vcf_s3 = join(work_dir, 'samples/a/b/variants/SNVs/snvs_NC_045512.2.snpSIFT.txt')

        run(f'sh -c "sed "s/MN908947.3/NC_045512.2/g" {vcf_s0} > {vcf_s1}"', params)
        # run('java -jar /opt/snpEff/snpEff.jar download -v NC_045512.2', params)
        run(f'sh -c "java -Xmx8g -jar /opt/snpEff/snpEff.jar NC_045512.2 {vcf_s1} > {vcf_s2}"', params)
        run(f'sh -c "cat {vcf_s2} | /opt/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /opt/snpEff/SnpSift.jar ' + 
            f' extractFields - CHROM POS REF ALT AF DP "ANN[*].IMPACT" "ANN[*].FEATUREID" "ANN[*].EFFECT" ' + 
            f' "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].AA_POS" "ANN[*].GENE" > {vcf_s3}"', params)

        smart_copy(join(work_dir, 'samples/a/b/alignments'), join(output_dir, 'alignments'))
        smart_copy(join(work_dir, 'samples/a/b/variants'), join(output_dir, 'variants'))
        smart_copy(join(work_dir, 'samples/a/b/extracted_data'), join(output_dir, 'extracted_data'), exclude=['*.gz', '*.fasta'])
        smart_copy(join(work_dir, 'samples/a/b/preprocessed_data'), join(output_dir, 'preprocessed_data'), exclude=['*.gz', '*.fasta'])
        smart_copy(join(work_dir, 'samples/a/b/raw_data'), join(output_dir, 'raw_data'), exclude=['*.gz', '*.fasta'])
        smart_copy(join(work_dir, 'samples/a/b/references'), join(output_dir, 'references'), exclude=['*.gz', '*.fasta'])
    
    # Assign parallel workers
    input_dirs = glob(join(args.input_dirs, '*'))
    if args.single_subject != '':
        input_dirs = [join(args.input_dirs, args.single_subject)]
        print(f"WARNING: only running --single_subject {args.single_subject}")
    results =  []
    for input_dir in input_dirs:
        output_dir = join(args.output_dirs, basename(input_dir))
        if exists(output_dir):
            if not args.force:
                print(f"WARNING: skipping {output_dir} since it already exists. Overwrite with --force.")
                continue
            else:
                smart_remove(output_dir)
        results.append(run_worker(abspath(input_dir), abspath(output_dir), params))

    for r in results:
        r.result()
