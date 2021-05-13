# Massively Parallel and Portable Genomic Sequence Analysis (Mappgene)

Mappgene is a genomic sequence analysis workflow designed for high-performance computing. It currently wraps V-pipe (https://github.com/cbg-ethz/V-pipe) with a collection of useful scripts for deployment in almost any Linux environment. 

## Quick Install

```
pip3 install --user parsl
git lfs install
git clone https://github.com/LLNL/mappgene.git
```

## Running

### Organizing your data

1. Copy your reference genome in `fasta` format to `vpipe_files/references/`.
2. Inspect and configure `vpipe_files/vpipe.config` for your run (e.g.,
   `reference = `, `trim_percent_cutoff = `, `threads = `)
3. Organize your input files (typically gzip compressed `fastq` formatted
   paired-end reads). `mappgene.py` expects the following layout, where sample
   names are taken from the subdirectories (i.e., foo, bar, baz):

    ```
    /path/to/input_dirs
    |-- foo
    |   |-- foo_R1.fastq.gz
    |   `-- foo_R2.fastq.gz
    |-- bar
    |   |-- bar_R1.fastq.gz
    |   `-- bar_R2.fastq.gz
    |-- baz
    |   |...
    ...
    ```

### Run

<b>Specify arguments</b>

Example:

```bash
python3 mappgene.py \
    --input_dirs /path/to/input_dirs \
    --output_dirs /path/to/output_dirs \
    --read_length 130 \
    --nnodes 16 \
    --walltime 12:59:59
```

**OR**

`python3 mappgene.py <config_json>`

<b>More info</b>

`python3 mappgene.py --help`

## How it works

1. V-pipe and other software have been pre-installed in a container (https://github.com/hpcng/singularity).
2. User parameters are parsed from command-line or configuration JSON (e.g. `configs/example/catalyst.json`). Note that every parameter has a default value (see `mappgene.py`).
3. Software tasks are then distributed to compute nodes in parallel (https://github.com/Parsl/parsl).

## License

Mappgene is distributed under the terms of the BSD-3 License.

LLNL-CODE-821512
