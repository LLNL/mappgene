#!/usr/bin/env python3
import argparse,parsl,os,sys,glob,shutil
from subscripts import *

parser = argparse.ArgumentParser()

parser.add_argument('--inputs', '-i', nargs='+', default='inputs/*.fastq.gz',
    help='Paths to FASTQ input file(s).')

parser.add_argument('--outputs', '-o', default='outputs/',
    help='Path to output directory.')

parser.add_argument('--read_length', default=250,
    help='Read length in sample.tsv (see cbg-ethz.github.io/V-pipe/tutorial/sars-cov2).')

parser.add_argument('--nnodes', '-n', default=1,
    help='Number of nodes.')

parser.add_argument('--bank', '-b', default='asccasc',
    help='Bank to charge for jobs.')

parser.add_argument('--partition', '-p', default='pbatch',
    help='Scheduler partition to assign jobs.')

parser.add_argument('--container', default='container/image.sif',
    help='Path to Singularity container image.')

parser.add_argument('--walltime', '-t', default='11:59:00',
    help='Walltime in format HH:MM:SS.')

parser.add_argument('--local', '-l', action='store_true',
    help='Run only on local threads (for debugging).')

parser.add_argument('--flux', action='store_true',
    help='Use the Flux scheduler instead of Parsl\'s native scheduler.')

parser.add_argument('--no_ivar', action='store_true',
    help='Disable ivar step.')

parser.add_argument('--no_vpipe', action='store_true',
    help='Disable vpipe step.')

parser.add_argument('--no_overwrite', action='store_true',
    help='Don\'t overwrite existing outputs.')

args = parser.parse_args()

if __name__ == '__main__':

    if isinstance(args.inputs, str):
        args.inputs = glob(args.inputs)

    deinterleaved_inputs = []
    merged_inputs = []
    for read1 in args.inputs:
        if '_R1.' in read1:
            read2 = read1.replace('_R1.', '_R2.')
            deinterleaved_inputs.append((read1, read2))
        elif '_R2.' in read1:
            pass
        else:
            merged_inputs.append(read1)

    # Copy V-pipe repo as main working directory
    tmp_dir = os.path.abspath('tmp')
    git_dir = os.path.join(tmp_dir, 'vpipe')
    smart_remove(tmp_dir)
    smart_mkdir(tmp_dir)
    git_params = {'sdir':git_dir, 'container':args.container}
    run(f'cp -rf /opt/vpipe {git_dir}', git_params)
    smart_copy('extra_files/', git_dir)
    run(f'cd {git_dir} && sh init_project.sh || true', git_params)

    params = {
        'container': os.path.abspath(args.container),
        'git_dir': os.path.abspath(git_dir),
        'read_length': args.read_length,
        'script_dir': os.path.dirname(os.path.abspath(__file__)),
    }

    if args.local:
        executor = parsl.executors.ThreadPoolExecutor(label="worker")
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

    config = parsl.config.Config(executors=[executor])
    parsl.set_stream_logger()
    parsl.load(config)

    results =  []
    for merged_input in merged_inputs:
        read1 = replace_extension(merged_input, '_R1.fastq.gz')
        read2 = replace_extension(merged_input, '_R2.fastq.gz')
        print(read1, read2)
        results.append(run_deinterleave(merged_input))
    for r in results:
        r.result()



    # results =  []
    # for input_pair in deinterleaved_inputs:
    #     output_dir = join(args.output_dirs, k)
    #     if exists(output_dir):
    #         if not args.force:
    #             print(f"WARNING: skipping {output_dir} since it already exists. Overwrite with --force.")
    #             continue
    #         else:
    #             smart_remove(output_dir)
    #     results.append(run_vpipe(inputs_dict[k], os.path.abspath(output_dir), params))

    # for r in results:
    #     r.result()
