#!/usr/bin/env python3
import argparse,parsl,os,sys,glob,shutil
from os.path import *
from subscripts import *

parser = argparse.ArgumentParser()

parser.add_argument('--inputs', '-i', nargs='+', default='inputs/*.fastq.gz',
    help='Paths to FASTQ input file(s).')

parser.add_argument('--outputs', '-o', default='outputs/',
    help='Path to output directory.')

parser.add_argument('--read_length', default=120,
    help='Read length in sample.tsv (see cbg-ethz.github.io/V-pipe/tutorial/sars-cov2).')

parser.add_argument('--container', default='container/image.sif',
    help='Path to Singularity container image.')

parser.add_argument('--walltime', '-t', default='11:59:00',
    help='Walltime in format HH:MM:SS.')

parser.add_argument('--no_ivar', action='store_true',
    help='Disable ivar step.')

parser.add_argument('--no_vpipe', action='store_true',
    help='Disable vpipe step.')

scheduler_group = parser.add_mutually_exclusive_group()

scheduler_group.add_argument('--slurm', '-s', action='store_true',
    help='Use the Slurm scheduler.')

scheduler_group.add_argument('--flux', action='store_true',
    help='Use the Flux scheduler.')

parser.add_argument('--nnodes', '-n', default=1,
    help='Slurm/Flux: number of nodes.')

parser.add_argument('--bank', '-b', default='asccasc',
    help='Slurm/Flux: bank to charge for jobs.')

parser.add_argument('--partition', '-p', default='pbatch',
    help='Slurm/Flux: partition to assign jobs.')

args = parser.parse_args()

if __name__ == '__main__':

    # Copy V-pipe repo as main working directory
    tmp_dir = abspath('tmp')
    vpipe_dir = join(tmp_dir, 'vpipe')
    base_params = {
        'container': abspath(args.container),
        'work_dir': tmp_dir,
        'read_length': args.read_length,
    }
    smart_remove(tmp_dir)
    smart_mkdir(tmp_dir)
    
    run(f'cp -rf /opt/vpipe {vpipe_dir}', base_params)
    smart_copy('extra_files', vpipe_dir)
    run(f'cd {vpipe_dir} && sh init_project.sh || true', base_params)
    update_permissions(base_params)

    if isinstance(args.inputs, str):
        args.inputs = glob(args.inputs)
    
    all_params = {}

    # Copy reads to subject directory
    for f in args.inputs:
        out1 = f.replace('_R2.', '_R1.')
        out2 = f.replace('_R1.', '_R2.')
        if '_R1.' in f and out2 not in args.inputs:
            raise Exception(f'Missing paired read: {out2}')
        if '_R2.' in f and out1 not in args.inputs:
            raise Exception(f'Missing paired read: {out1}')

        subject = replace_extension(basename(f).replace('_R1.', '.').replace('_R2.', '.'))
        subject_dir = abspath(join(args.outputs, subject))
        raw_dir = join(vpipe_dir, 'samples/a/b/raw_data')
        read = join(raw_dir, basename(f))

        if not subject in all_params:
            smart_remove(raw_dir)
            params = base_params.copy()
            params['work_dir'] = subject_dir
            params['stdout'] = join(subject_dir, 'mappgene.stdout')
            all_params[subject] = params

        smart_copy(tmp_dir, subject_dir)
        smart_copy(f, read)

    if args.slurm:
        executor = parsl.executors.HighThroughputExecutor(
            label="worker",
            address=parsl.addresses.address_by_hostname(),
            provider=parsl.providers.SlurmProvider(
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
    elif args.flux:
        executor = parsl.executors.FluxExecutor(
            label="worker",
            flux_path="/usr/global/tools/flux/toss_3_x86_64_ib/flux-c0.28.0.pre-s0.17.0.pre/bin/flux",
            provider=parsl.providers.SlurmProvider(
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
    else:
        executor = parsl.executors.ThreadPoolExecutor(label="worker")

    config = parsl.config.Config(executors=[executor])
    parsl.set_stream_logger()
    parsl.load(config)

    if not args.no_ivar:
        results =  []
        for params in all_params.values():
            results.append(run_ivar(params))
        for r in results:
            r.result()

    if not args.no_vpipe:
        results =  []
        for params in all_params.values():
            results.append(run_vpipe(params))
        for r in results:
            r.result()
