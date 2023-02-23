mappgene
===========
[![PyPI version](https://badge.fury.io/py/mappgene.svg)](https://badge.fury.io/py/mappgene) [![DOI](https://zenodo.org/badge/367142630.svg)](https://zenodo.org/badge/latestdoi/367142630)

mappgene is a SARS-CoV-2 genomic sequence analysis pipeline designed for high-performance parallel computing. It mainly wraps iVar (https://github.com/andersen-lab/ivar) and LoFreq (https://github.com/CSB5/lofreq) with a collection of useful scripts for deployment in almost any Linux environment.

Inputs: `.fastq.gz`

Outputs: `.vcf` and `.variants.tsv`


## Quick Setup

#### Requirements

* Python 3.7+
* [Singularity](https://sylabs.io/guides/3.5/user-guide/index.html)

#### Installation
```
pip3 install mappgene
singularity pull library://khyox/mappgene/image.sif:latest
```

#### Usage
```
mappgene <SUBJECT.FASTQ.GZ>
```

## Instructions

#### Example Testing
Check that mappgene works on your system by running the example input data, sourced from [here](https://github.com/cbg-ethz/V-pipe/tree/master/tests/data/sars-cov-2/pos_MN908947_3_1/20200729/raw_data).
```
mappgene --test
```

#### Multiple subjects
You can specify multiple subjects with specific paths or Unix-style globbing
```
mappgene <SUBJECT1.FASTQ.GZ> <SUBJECT2.FASTQ.GZ> <SUBJECT3.FASTQ.GZ>
mappgene <SUBJECT_DIR>/*.fastq.gz
```

#### Deinterleaved subjects
If there are two subjects with matching names that end in `_R1.fastq.gz` and `_R2.fastq.gz`, mappgene will assume they are a deinterleaved pair.
```
mappgene <SUBJECT>_R1.fastq.gz <SUBJECT>_R1.fastq.gz
```

#### Slurm scheduling
Multiple subjects can be run on distributed systems using Slurm or Flux.
```
mappgene --slurm -n 1 -b mybank -p mypartition <SUBJECT.FASTQ.GZ>
```

#### Additional options
```
mappgene --help
```

License
-------
Mappgene is distributed under the terms of the BSD-3 License.

LLNL-CODE-821512

____
You may be interested in [MappgeneSummary](https://github.com/LLNL/mappgenesummary), a package for the analysis and summarization of mappgene's results.

If you use mappgene in your research, please cite the paper. Thanks!

Kimbrel J, Moon J, Avila-Herrera A, Martí JM, Thissen J, Mulakken N, Sandholtz SH, Ferrell T, Daum C, Hall S, Segelke B, Arrildt KT, Messenger S, Wadford DA, Jaing C, Allen JE, Borucki MK. Multiple Mutations Associated with Emergent Variants Can Be Detected as Low-Frequency Mutations in Early SARS-CoV-2 Pandemic Clinical Samples. Viruses. 2022; 14(12):2775. https://doi.org/10.3390/v14122775
____
