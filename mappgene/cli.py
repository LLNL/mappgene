#!/usr/bin/env python3
import argparse,parsl,os,sys,glob,shutil
from os.path import *
from mappgene.subscripts import *

script_dir = abspath(os.path.dirname(os.path.realpath(__file__)))
cwd = abspath(os.getcwd())

def parse_args(args):


    parser = argparse.ArgumentParser()

    if not '--test' in sys.argv:
        parser.add_argument('inputs', nargs='+',
            help='Paths to FASTQ input file(s).')

    parser.add_argument('--test', action='store_true',
        help='Test using the example inputs.')

    workflow_group = parser.add_mutually_exclusive_group(required=True)

    workflow_group.add_argument('--ivar', action='store_true',
        help='Run the iVar workflow.')

    workflow_group.add_argument('--vpipe', action='store_true',
        help='Run the V-pipe workflow.')

    parser.add_argument('--outputs', '-o', default='mappgene_outputs/',
        help='Path to output directory.')

    parser.add_argument('--container', default=join(cwd, 'image.sif'),
        help='Path to Singularity container image.')

    parser.add_argument('--read-length', default=130,
        help='V-pipe: read length in sample.tsv (see cbg-ethz.github.io/V-pipe/tutorial/sars-cov2).')

    parser.add_argument('--variant_frequency', default=0.01,
        help='iVar: variant frequency cutoff.')

    scheduler_group = parser.add_mutually_exclusive_group()

    scheduler_group.add_argument('--slurm', action='store_true',
        help='Use the Slurm scheduler.')

    scheduler_group.add_argument('--flux', action='store_true',
        help='Use the Flux scheduler.')

    parser.add_argument('--nnodes', '-n', default=1,
        help='Slurm/Flux: number of nodes.')

    parser.add_argument('--bank', '-b', default='asccasc',
        help='Slurm/Flux: bank to charge for jobs.')

    parser.add_argument('--partition', '-p', default='pbatch',
        help='Slurm/Flux: partition to assign jobs.')

    parser.add_argument('--walltime', '-t', default='11:59:00',
        help='Slurm/Flux: walltime in format HH:MM:SS.')

    return parser.parse_args()

def main():

    args = parse_args(sys.argv[1:])

    # Copy V-pipe repo as main working directory
    tmp_dir = join(cwd, 'tmp')
    vpipe_dir = join(tmp_dir, 'vpipe')
    base_params = {
        'container': abspath(args.container),
        'work_dir': tmp_dir,
        'read_length': args.read_length,
        'variant_frequency': args.variant_frequency,
        'stdout': abspath(join(args.outputs, 'mappgene.stdout')),
    }

    if shutil.which('singularity') is None:
        raise Exception(f"Missing Singularity executable in PATH.\n\n" +
            f"Please ensure Singularity is installed: https://sylabs.io/guides/3.0/user-guide/installation.html")

    if not exists(base_params['container']):
        raise Exception(f"Missing container image at {base_params['container']}\n\n" +
            f"Either specify another image with --container\n\n" +
            f"Or build the container with the recipe at: {join(script_dir, 'data/container/recipe.def')}\n\n" +
            f"Or download the container at: https://www.dropbox.com/s/ymsfn9z7v3utqe0/image.sif?dl=1\n")

    smart_remove(tmp_dir)
    smart_mkdir(tmp_dir)
    
    run(f'cp -rf /opt/vpipe {vpipe_dir}', base_params)
    smart_copy(join(script_dir, 'data/extra_files'), tmp_dir)

    run(f'cd {vpipe_dir} && sh init_project.sh || true', base_params)
    update_permissions(tmp_dir, base_params)

    if args.test:
        args.inputs = join(script_dir, 'data/example_inputs/*.fastq.gz')
    
    if isinstance(args.inputs, str):
        args.inputs = glob(args.inputs)

    all_params = {}

    # Copy reads to subject directory
    for input_read in args.inputs:
        pair1 = input_read.replace('_R2.', '_R1.')
        pair2 = input_read.replace('_R1.', '_R2.')
        if '_R1.' in input_read and pair2 not in args.inputs:
            raise Exception(f'Missing paired read: {pair2}')
        if '_R2.' in input_read and pair1 not in args.inputs:
            raise Exception(f'Missing paired read: {pair1}')

        subject = replace_extension(basename(input_read).replace('_R1.', '.').replace('_R2.', '.'))
        subject_dir = abspath(join(args.outputs, subject))

        if not subject in all_params:
            smart_copy(tmp_dir, subject_dir)
            params = base_params.copy()
            params['work_dir'] = subject_dir
            params['input_reads'] = [input_read]
            params['stdout'] = join(subject_dir, 'worker.stdout')
            all_params[subject] = params
        else:
            all_params[subject]['input_reads'].append(input_read)

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

    if args.ivar:
        results =  []
        for params in all_params.values():
            results.append(run_ivar(params))
        for r in results:
            r.result()

    elif args.vpipe:
        results =  []
        for params in all_params.values():
            results.append(run_vpipe(params))
        for r in results:
            r.result()

if __name__ == '__main__':
    main()
