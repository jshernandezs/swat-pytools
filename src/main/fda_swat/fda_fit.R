#!/usr/bin/env Rscript
#### import packages ####
suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c('-i','--input'), action = 'store', default = './Data/observed.csv', type = 'character',
              help = 'Directory of.csv file with time series (including file name) [default %default]',
              metavar = ""),
  make_option(c('-o','--output'), action = 'store', default = './Data/observed_fmodel.rds', type = 'character',
              help = 'Directory of.rds file with Fourier series model (including file name) [default %default]',
              metavar = ""),
  make_option(c('-n','--nbasis'), action = 'store', default = 365, type = 'integer',
              help = 'Number of basis for fitting models [default %default]',
              metavar = ""),
  make_option(c('-l','--log10lambda'), action = 'store', default = 'seq(-4,2,0.25)', type = 'character',
              help = 'R expression to build the numeric array of log10 lambda penalizing coefficients to evaluate for cross validation [default %default]',
              metavar = ""),
  make_option(c('-p','--plot'), action = 'store_true', default = FALSE,
              help = 'Do plot results [default %default]',
              metavar = ""),
  make_option(c('-f','--fast'), action = 'store_true', default = FALSE,
              help = 'Do run in fast mode [default %default] Warning: only works for complete/regular time series',
              metavar = "")
  )

opt = parse_args(OptionParser(option_list = option_list))

#### input parameters ####
obs.file    <- opt$input  # csv file with time series
K           <- opt$nbasis # number of basis
loglam      <- eval(parse(text = opt$log10lambda)) # sequence of logaritmic lambda coefs. to evaluate
opt_plot    <- opt$plot   # boolean variable for plotting results
output.file <- opt$output # rds output file with model fitting results
opt.fast    <- opt$fast   # boolean for running in fast mode

#### main script ####
# opt.fast <- TRUE
source('./fda_utilities.R')
if (opt.fast){
  run_fit_fast(obs.file,K,loglam,opt_plot,output.file)
} else {
  run_fit(obs.file,K,loglam,opt_plot,output.file)
}