#!/usr/bin/env Rscript
#### import packages ####
suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c('-s','--simulated'), action = 'store', default = './Data/simulated.csv', type = 'character',
              help = 'Directory of.csv file with simulated time series (including file name) [default %default]',
              metavar = ""),
  make_option(c('-o','--observed'), action = 'store', default = './Data/observed_fmodel.rds', type = 'character',
              help = 'Directory of.rds file with Fourier series model (including file name) [default %default]',
              metavar = ""),
  make_option(c('-p','--plot'), action = 'store_true', default = FALSE,
              help = 'Do plot results [default %default]',
              metavar = ""),
  make_option(c('-z','--zero'), action = 'store_true', default = FALSE,
              help = 'Do set lambda to zero [default %default] (Note: only works for fast mode)',
              metavar = ""),
  make_option(c('-f','--fast'), action = 'store_true', default = FALSE,
              help = 'Do run in fast mode [default %default]',
              metavar = ""),
  make_option(c('-n','--new'), action = 'store_true', default = FALSE,
              help = 'Do fit new lambda for simulated time series [default %default]',
              metavar = ""),
  make_option(c('-l','--log10lambda'), action = 'store', default = 'seq(-4,2,0.25)', type = 'character',
              help = 'R expression to build the numeric array of log10 lambda penalizing coefficients to evaluate for cross validation [default %default]',
              metavar = "")
)

opt = parse_args(OptionParser(option_list = option_list))

#### input parameters ####
obs.saved.file <- opt$obs
sim.file       <- opt$sim
opt.plot       <- opt$plot
opt.zero       <- opt$zero
opt.fast       <- opt$fast
opt.new        <- opt$new
loglam         <- eval(parse(text = opt$log10lambda))

#### main script ####
source('./fda_utilities.R')
of <- get_metric(obs.saved.file, sim.file, opt.plot, opt.zero, opt.fast, opt.new, loglam)
cat(sprintf('%.12f',of))