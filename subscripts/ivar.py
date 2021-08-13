#!/usr/bin/env python3
from parsl.app.app import python_app

@python_app(executors=['worker'], cache=True)
def run_ivar(params, inputs=[]):
    pass