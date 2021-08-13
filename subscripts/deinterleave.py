#!/usr/bin/env python3
from parsl.app.app import python_app

@python_app(executors=['worker'], cache=True)
def run_deinterleave(merged_input, inputs=[]):
    read1 = merged_input.split('.')[0] + '_R1.fastq.gz'
    read2 = merged_input.split('.')[0] + '_R2.fastq.gz'
    pass