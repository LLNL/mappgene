mappgene
========

[![DOI](https://zenodo.org/badge/367142630.svg)](https://zenodo.org/badge/latestdoi/367142630)

mappgene is a SARS-CoV-2 variant calling pipeline designed for high-performance
parallel computing. It mainly wraps iVar (https://github.com/andersen-lab/ivar)
and LoFreq (https://github.com/CSB5/lofreq) variant callers, plus
snpEff/snpSift annotation tools.

Inputs: short-read paired-end Illumina RTA3 sequencing data in fastq format
(gzip compressed) (e.g., `SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz`)

Outputs: variants in `.vcf` variant call format files, and snpEff/snpSift
tabular text output

## Quick Start

```bash
singularity pull library://khyox/mappgene/image.sif:latest
git clone https://github.com/LLNL/mappgene.git
pip install git+file:///absolute/path/to/mappgene
mappgene --container image.sif --outputs outputs samples/*fastq.gz
```

## Requirements

* Python 3.7+
* [Singularity](https://sylabs.io/guides/3.5/user-guide/index.html)

## Installation

We recommend installing mappgene to a python3 virtualenv. We have found
it useful to install in "editable" mode to easily customize and modify
mappgene.

```bash
python3 -m venv mg
source mg/bin/activate
pip install -e git+file:///absolute/path/to/mappgene
```

Don't forget to download the corresponding singularity container!
It contains the pipeline components and dependencies.

- iVar v1.3.2
- LoFreq 2.1.5
- See `mappgene/data/container/recipe.def` for more details

```bash
singularity pull library://khyox/mappgene/image.sif:latest
```

Or go to <https://cloud.sylabs.io/library/khyox/mappgene/image.sif> and click
"Download".

### Example Testing

Check that mappgene works on your system by running the example input data,
sourced from
[here](https://github.com/cbg-ethz/V-pipe/tree/master/tests/data/sars-cov-2/pos_MN908947_3_1/20200729/raw_data).

```bash
mappgene --test
```

## Usage

See `mappgene -h` for a list of options and detailed usage. In short:

```
mappgene [OPTIONS] <SAMPLE1_R1.fastq.gz> <SAMPLE1_R2.fastq.gz> [SAMPLEN_R1.fastq.gz ...]
```

Key options:

- `--container`: tell mappgene where the singularity container is
- `--slurm`: use the slurm scheduler for processing samples in parallel
- `--use_full_node`: 1 sample per node
- `--primers_bp`: specify a bundled primer set to use
- `--depth_cap`: sets `lofreq call -d` value (read no more than this many reads per position)
- `--read_cutoff_bp`: sets `ivar trim -m` value (remove reads smaller than this after trimming)
- `--variant_frequency`: sets `ivar variants -t` value (do not call variants below this frequency)

## Instructions

### Process multiple samples in parallel

You can specify multiple samples with specific paths or Unix-style globbing.
Reads must be gzip compressed (`.fastq.gz`).

If there are two input filenames with a matching sample name, plus `_R1`
and `_R2`, then mappgene will assume they are a deinterleaved pair. For
deinterleaved reads, you must ensure:

- first and second read files contain `_R1` and `_R2`, respectively
- there is only one pair of read files per sample (no orphans, no multi-run
  samples, samples from multi-sample subjects are treated separately)
- read pairs appear together on the command line when expanding shell globs

```
mappgene <SAMPLE1>_R1.fastq.gz <SAMPLE1>_R2.fastq.gz [<SAMPLE2>_R1.fastq.gz <SAMPLE2>_R2.fastq.gz ...]
```

#### Interleaved samples

If the input filenames do not contain `_R1` or `_R2`, mappgene will probably
interpret the inputs as interleaved samples, and automatically deinterleave
them during processing.

```
mappgene <SAMPLE1.FASTQ.GZ> <SAMPLE2.FASTQ.GZ> <SAMPLE3.FASTQ.GZ>
mappgene <SAMPLE_DIR>/*.fastq.gz
```

#### Slurm scheduling

Multiple subjects can be run in parallel on HPC systems using the Slurm job
scheduler.

```
mappgene --slurm -n 1 -b mybank -p mypartition <SAMPLE.FASTQ.GZ>
```

### Output

By default, results will be in `mappgene_outputs/<SAMPLE>`, or
wherever specified by `--outputs`.

Key output files:

```bash
mappgene_outputs/
  <SAMPLE>/
    worker.stdout  # (log file capturing stdout)
    ivar_outputs/
      <SAMPLE>.ivar.snpEff.vcf  # (ivar variant calls)
      <SAMPLE>.ivar.snpSift.txt
      <SAMPLE>.ivar.lofreq.snpEff.vcf  # (lofreq variant calls)
      <SAMPLE>.ivar.lofreq.snpSift.txt
```

## Known bugs and quirks

- The `--flux` option to use the Flux scheduler is broken.
- Only `fastp` respects mappgene's `--threads` option. `bwa mem` and `lofreq`
  use different numbers of threads. Additionally, if `--use_full_node` is not
  specified, mappgene will try to run multiple samples per node.
- `fastp`, `snpEff`, `snpSift` write to disk outside of the current working
  directory—by default to the user's home—which may be on a different
  filesystem not intended for parallel I/O.  This output is typically not used,
  but it can still clobber or corrupt existing files, or impact cluster
  performance for all users.
- Viral-recon's `ivar_variants_to_vcf.py` attempts to group consecutive SNPs in
  the same codon into single multinucleotide variants. Previous combinations of
  `ivar` and viral-recon script versions have introduced runtime errors
  (resolved by using patched `dev` branch versions).
- Running multiple instances of mappgene in the same working directory is **not
  recommended**, as a common temporary directory is used, resulting in scary
  warnings, possible errors, and potentially corrupted runs.
- `mappgene` assumes quality scores are RTA3 "score category" labels (error,
  low, medium, high) for base calls rather than a continuous numeric score.
  Although the quality score is supposed to reflect an average score for each
  base call category, mappgene takes a conservative approach and adjusts medium
  and high scores to the lower bound of those categories, i.e., Q37->Q30 and
  Q25->Q20. This affects `lofreq`, a quality-aware variant caller, by making it
  require more evidence (i.e., depth of read coverage) to call variants vs.
  sequencing errors.

License
-------

Mappgene is distributed under the terms of the BSD-3 License.

LLNL-CODE-821512

____
You may be interested in
[MappgeneSummary](https://github.com/LLNL/mappgenesummary), a package for the
analysis and summarization of mappgene's results.

If you use mappgene in your research, please cite the paper. Thanks!

Kimbrel J, Moon J, Avila-Herrera A, Martí JM, Thissen J, Mulakken N, Sandholtz
SH, Ferrell T, Daum C, Hall S, Segelke B, Arrildt KT, Messenger S, Wadford DA,
Jaing C, Allen JE, Borucki MK. Multiple Mutations Associated with Emergent
Variants Can Be Detected as Low-Frequency Mutations in Early SARS-CoV-2
Pandemic Clinical Samples. Viruses. 2022; 14(12):2775.
https://doi.org/10.3390/v14122775
____
