#!/usr/bin/env python3
"""
metrics from functional data analysis

"""
import subprocess as sp
import os


class MetricFda:
    def __init__(self, fda_dir=None):
        if fda_dir is None:
            self.fda_dir = os.path.abspath('./fda_swat')
        else:
            self.fda_dir = os.path.abspath(fda_dir) 
        
        self.metric = []

    def fda_fit(self, file, output_file, k=365, loglam='seq(-4,2,0.25)', opt_plot=False, opt_fast=False):
        line = ['./fda_fit.R', '-i', file, '-o', output_file, '-n', str(k), '-l', loglam]
        if opt_plot:
            line.append('-p')
        if opt_fast:
            line.append('-f')
        sp.run(line, cwd=self.fda_dir, check=True)

    def fda_metric(self, file, output_file, opt_plot=False, opt_zero=False, opt_fast=False, opt_new=False,
                   loglam='seq(-4,2,0.25)'):
        line = ['./fda_metric.R', '-s', file, '-o', output_file, '-l', loglam]
        if opt_plot:
            line.append('-p')
        if opt_fast:
            line.append('-f')
        if opt_zero:
            line.append('-z')
        if opt_new:
            line.append('-n')
        p = sp.run(line, cwd=self.fda_dir, check=True, stdout=sp.PIPE)
        metric = float(p.stdout)
        self.metric = metric
        return metric
