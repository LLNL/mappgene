mappgene
===========
[![PyPI version](https://badge.fury.io/py/mappgene.svg)](https://badge.fury.io/py/mappgene) [![DOI](https://zenodo.org/badge/367142630.svg)](https://zenodo.org/badge/latestdoi/367142630)

mappgene is a SARS-CoV-2 genomic sequence analysis pipeline designed for high-performance parallel computing. It wraps V-pipe (https://github.com/cbg-ethz/V-pipe) and iVar (https://github.com/andersen-lab/ivar) with a collection of useful scripts for deployment in almost any Linux environment.

Inputs: `.fastq.gz`

Outputs: `.vcf` and `.snpSIFT.txt`


## Quick Setup

#### Requirements

* Python 3.7+
* [Singularity](https://sylabs.io/guides/3.5/user-guide/index.html)

#### Installation
```
pip3 install mappgene
singularity pull image.sif library://avilaherrera1/mappgene/image.sif:latest
```

#### Usage
```
mappgene --ivar <SUBJECT.FASTQ.GZ>
```

## Instructions

#### Example Testing
Check that mappgene works on your system by running the example input data, sourced from [here](https://github.com/cbg-ethz/V-pipe/tree/master/tests/data/sars-cov-2/pos_MN908947_3_1/20200729/raw_data).
```
mappgene --ivar --test
```

#### Multiple subjects
You can specify multiple subjects with specific paths or Unix-style globbing
```
mappgene --ivar <SUBJECT1.FASTQ.GZ> <SUBJECT2.FASTQ.GZ> <SUBJECT3.FASTQ.GZ>
mappgene --ivar <SUBJECT_DIR>/*.fastq.gz
```

#### Deinterleaved subjects
If there are two subjects with matching names that end in `_R1.fastq.gz` and `_R2.fastq.gz`, mappgene will assume they are a deinterleaved pair.
```
mappgene --ivar <SUBJECT>_R1.fastq.gz <SUBJECT>_R1.fastq.gz
```

#### Run V-pipe
You can run either iVar or V-pipe.
```
mappgene --vpipe <SUBJECT.FASTQ.GZ>
```

#### Slurm scheduling
Multiple subjects can be run on distributed systems using Slurm or Flux.
```
mappgene --ivar --slurm -n 1 -b mybank -p mypartition <SUBJECT.FASTQ.GZ>
```

#### Additional options
```
mappgene --help
```

License
-------
Mappgene is distributed under the terms of the BSD-3 License.

LLNL-CODE-821512
