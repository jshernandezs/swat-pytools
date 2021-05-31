#### import packages ####
suppressPackageStartupMessages(require(fda))
suppressPackageStartupMessages(require(dplyr))

#### set of functions for performing Funtional Data Analysis (FDA) for regular time series ####
# determine days in a year
days_in_year  <- function(year) {
  365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
}
# obtain finite Fourier series for observed data using generalized cross-validation
best_fit_fast <- function(ts, loglam, K, fdnames){
  dayrange <- c(0, 365)
  # functional data object settings
  flowbasis <- create.fourier.basis(dayrange, K)
  Lcoef <- c(0,(2*pi/diff(dayrange))^2,0)
  harmaccelLfd <- vec2Lfd(Lcoef, dayrange)
  nlam <- length(loglam)
  gcvsave <- rep(NA,nlam)
  modelsave <- list()
  pb <- txtProgressBar(min = 1, max =nlam, style = 3)
  for (ilam in 1:nlam){
    # set smoothing parameters
    lambda  <- 10^loglam[ilam]
    flowPar <- fdPar(flowbasis, harmaccelLfd, lambda)
    smoothlist <- smooth.basis(day.5,ts,flowPar, fdnames = fdnames)
    # write-out results
    gcvsave[ilam]     <- sum(smoothlist$gcv)
    modelsave[[ilam]] <- smoothlist
    # update progress bar
    setTxtProgressBar(pb, ilam)
  }
  close(pb)
  best_gcv    <- min(gcvsave)
  best_model  <- modelsave[[which(gcvsave == best_gcv[1])[1]]]
  best_lambda <- 10^(loglam[which(gcvsave == best_gcv[1])[1]])
  return(list(best_model = best_model, gcvsave = gcvsave, lambda = best_lambda))
}
# obtain models using a global lambda
fda_fit_fast  <- function(obs.file, K, loglam){
  # read time series
  data <- read.csv(obs.file, header = FALSE)
  # set date format
  data[,1] <- as.POSIXct(data[,1], format = '%Y-%m-%d')
  # remove additional day from leap years
  data <- data[format(data[,1], '%m-%d') != '02-29',]
  # get list of years in time series
  years <- unique(format(data[,1], '%Y'))
  # obtain best finite Fourier series for each year
  cat('Obtaining best Fourier models for designated time series...\n')
  ts <- numeric()
  for (year in years){
    ts.annual  <- data[format(data[,1],'%Y') %in% year, ] # annual time-series
    ts <- cbind(ts,ts.annual[,2])
  }
  fdnames    <- list('Time','Years' = years, 'Variable')
  best_list  <- best_fit_fast(ts,loglam,K,fdnames)
  return(best_list)
}
# plot best results
plot_fit_fast <- function(best_list, loglam, opt_plot){
  if (opt_plot){
    best_model <- best_list$best_model
    gcvsave    <- best_list$gcvsave
    # best model plot
    ypred <- eval.fd(day.5,best_model$fd)
    yobs  <- best_model$y
    for (i in 1:dim(yobs)[2]){
      par()
      plot(dayOfYear,yobs[,i], xlab = 'Date', ylab = 'Variable', pch = 20)
      lines(dayOfYear,ypred[,i])
    }
    par()
    plot(loglam,gcvsave, ylab = 'GCV criterion', xlab = 'log10(lambda)', type = 'b')
    # plot overall performace
    ypred <- as.vector(ypred)
    yobs  <- as.vector(yobs)
    par(mfrow = c(1,1), pty = 's')
    plot(ypred,yobs,xlab = 'Predicted variable', ylab = 'Observed variable',
         xlim = range(ypred, yobs), ylim = range(ypred, yobs), pch = 20)
    abline(coef = c(0,1))
  }
}
# save results for fast fit
save_fit_fast <- function(best_list, output.file){
  # modeled and simulated outputs
  best_model <- best_list$best_model
  gcvsave    <- best_list$gcvsave
  lambda     <- best_list$lambda
  ypred <- eval.fd(day.5,best_model$fd)
  yobs  <- best_model$y
  # save output for each year
  years <- best_model$fd$fdnames$Years
  colnames(ypred) <- colnames(yobs) <- years
  # obtain overall performace
  ypred2 <- as.vector(ypred)
  yobs2  <- as.vector(yobs)
  df   <- data.frame(yobs = yobs2, ypred = ypred2)
  beta <- lm(yobs ~ ypred, df)$coefficients
  R2   <- cor(yobs2,ypred2)^2
  if (beta[2] <= 1){
    bR2 <- beta[2]*R2
  } else {
    bR2 <- R2/beta[2]
  }
  cat("\n")
  cat(sprintf('Overall R2 is %.4f',R2))
  cat("\n")
  cat(sprintf('Overall bR2 is %.4f',bR2))
  cat("\n")
  cat(sprintf('Intercept is %.4f and slope is %.4f\n',beta[1],beta[2]))
  fmodel.set <- list(best_model = best_model$fd, gcvsave = gcvsave, lambda = lambda, ypred = ypred, yobs = yobs, performance = list(R2 = R2, bR2 = bR2, beta = beta))
  saveRDS(fmodel.set, file = output.file)
}
# fit Fourier models using a global lambda (Note: only works for complete/regular time series)
run_fit_fast  <- function(obs.file,K,loglam,opt_plot,output.file){
  # fit the model
  best_list <- fda_fit_fast(obs.file, K, loglam)
  # plot the model
  plot_fit_fast(best_list, loglam, opt_plot)
  # save the results
  save_fit_fast(best_list, output.file)
}

#### set of functions for performing Funtional Data Analysis (FDA) for irregular time series ####
# obtain finite Fourier series for observed data using generalized cross-validation
best_fit <- function(ts, days, loglam, K){
  dayrange <- c(0, 365)
  # functional data object settings
  flowbasis <- create.fourier.basis(dayrange, K)
  Lcoef <- c(0,(2*pi/diff(dayrange))^2,0)
  harmaccelLfd <- vec2Lfd(Lcoef, dayrange)
  nlam <- length(loglam)
  gcvsave <- rep(NA,nlam)
  modelsave <- list()
  pb <- txtProgressBar(min = 1, max =nlam, style = 3)
  for (ilam in 1:nlam){
    # set smoothing parameters
    lambda  <- 10^loglam[ilam]
    flowPar <- fdPar(flowbasis, harmaccelLfd, lambda)
    smoothlist <- smooth.basis(days,ts,flowPar)
    # write-out results
    gcvsave[ilam]     <- sum(smoothlist$gcv)
    modelsave[[ilam]] <- smoothlist
    # update progress bar
    setTxtProgressBar(pb, ilam)
  }
  close(pb)
  best_gcv    <- min(gcvsave)
  best_model  <- modelsave[[which(gcvsave == best_gcv[1])[1]]]
  best_lambda <- 10^(loglam[which(gcvsave == best_gcv[1])[1]])
  return(list(best_model = best_model, gcvsave = gcvsave, lambda = best_lambda))
}
# obtain finite Fourier series for observed data using generalized cross-validation
fda_fit  <- function(obs.file, K, loglam){
  # read time series
  data <- read.csv(obs.file, header = FALSE)
  # set date format
  data[,1] <- as.POSIXct(data[,1], format = '%Y-%m-%d')
  # remove additional day from leap years
  data <- data[format(data[,1], '%m-%d') != '02-29',]
  # get list of years in time series
  years <- unique(format(data[,1], '%Y'))
  # obtain best finite Fourier series for each year
  cat('Obtaining best Fourier models for designated time series...\n')
  best_lambda <- rep(NA,length(years))
  gcvsave <- days_list <- ts <- list()
  coefmat <- numeric()
  c <- 1
  for (year in years){
    ts.annual  <- data[format(data[,1],'%Y') %in% year, ] # annual time-series
    ndays      <- days_in_year(as.numeric(year))
    days       <- as.numeric(strftime(ts.annual[,1], format = "%j"))
    if (ndays == 366){
      days[days>59] <- days[days>59] - 1
    }
    days              <- days - 0.5
    best_list_i       <- best_fit(ts.annual[,2],days,loglam,K)
    best_model_i      <- best_list_i$best_model
    gcvsave[[year]]   <- best_list_i$gcvsave
    best_lambda[c]    <- best_list_i$lambda
    ts[[year]]        <- ts.annual
    days_list[[year]] <- days
    coefmat           <- cbind(coefmat,coef(best_model_i))
    c <- c + 1
  }
  fdnames    <- list('Time','Years' = years, 'Variable')
  flowbasis  <- best_model_i$fd$basis
  best_model <- fd(coefmat, flowbasis, fdnames)
  best_list  <- list(best_model = best_model, gcvsave = gcvsave, lambda = best_lambda, days_list = days_list, ts_obs = ts)
  return(best_list)
}
# plot best results
plot_fit <- function(best_list, loglam, opt_plot){
  if (opt_plot){
    best_model <- best_list$best_model
    gcvsave    <- best_list$gcvsave
    years      <- best_model$fdnames$Years
    ts         <- best_list$ts_obs
    days_list  <- best_list$days_list
    # best model plot
    ypred_list <- yobs_list <- list()
    for (i in 1:length(years)){
      tt     <- ts[[years[i]]][,1]
      tt.num <- days_list[[years[i]]]
      ypred  <- eval.fd(tt.num,best_model)[,i]
      yobs   <- ts[[years[i]]][,2]
      par(mfrow = c(1,2))
      plot(tt, yobs, xlab = 'Date', ylab = 'Variable', pch = 20)
      lines(tt, ypred)
      par()
      plot(loglam,gcvsave[[years[i]]], ylab = 'GCV criterion', xlab = 'log10(lambda)', type = 'b')
      ypred_list[[years[i]]] <- ypred
      yobs_list[[years[i]]]  <- yobs
    }
    # plot overall performace
    ypred <- unlist(ypred_list)
    yobs  <- unlist(yobs_list)
    par(mfrow = c(1,1), pty = 's')
    plot(ypred,yobs,xlab = 'Predicted variable', ylab = 'Observed variable',
         xlim = range(ypred, yobs), ylim = range(ypred, yobs), pch = 20)
    abline(coef = c(0,1))
  }
}
# save results for fast fit
save_fit <- function(best_list, output.file){
  # modeled and simulated outputs
  best_model <- best_list$best_model
  gcvsave    <- best_list$gcvsave
  lambda     <- best_list$lambda
  years      <- best_model$fdnames$Years
  ts_obs     <- best_list$ts_obs
  days_list  <- best_list$days_list
  ypred_list <- yobs_list <- list()
  for (i in 1:length(years)){
    year   <- years[i]
    tt.num <- days_list[[year]]
    ypred  <- eval.fd(tt.num,best_model)[,i]
    yobs   <- ts_obs[[year]][,2]
    ypred_list[[year]] <- ypred
    yobs_list[[year]]  <- yobs
  }
  # obtain overall performace
  ypred <- unlist(ypred_list)
  yobs  <- unlist(yobs_list)
  df    <- data.frame(yobs = yobs, ypred = ypred)
  beta  <- lm(yobs ~ ypred, df)$coefficients
  R2    <- cor(yobs,ypred)^2
  if (beta[2] <= 1){
    bR2 <- beta[2]*R2
  } else {
    bR2 <- R2/beta[2]
  }
  cat("\n")
  cat(sprintf('Overall R2 is %.4f',R2))
  cat("\n")
  cat(sprintf('Overall bR2 is %.4f',bR2))
  cat("\n")
  cat(sprintf('Intercept is %.4f and slope is %.4f\n',beta[1],beta[2]))
  save
  ypred <- ypred_list
  yobs  <- yobs_list
  fmodel.set <- list(best_model = best_model, gcvsave = gcvsave, lambda = lambda, 
                     ypred = ypred, yobs = yobs, days_list = days_list, performance = list(R2 = R2, bR2 = bR2, beta = beta))
  saveRDS(fmodel.set, file = output.file)
}
# fit Fourier models using a local lambda (Note: works for regular and irregular time series)
run_fit  <- function(obs.file,K,loglam,opt_plot,output.file){
  # fit the model
  best_list <- fda_fit(obs.file, K, loglam)
  # plot the model
  plot_fit(best_list, loglam, opt_plot)
  # save the results
  save_fit(best_list, output.file)
}

#### set of functions for computing fda-based performance metrics for model calibration ####
# obtain finite Fourier series for simulated data
model_fit  <- function(ts,days,lambda,flowbasis, fdnames = NULL){
  dayrange     <- flowbasis$rangeval
  Lcoef        <- c(0,(2*pi/diff(dayrange))^2,0)
  harmaccelLfd <- vec2Lfd(Lcoef, dayrange)
  flowPar <- fdPar(flowbasis, harmaccelLfd, lambda)
  model   <- smooth.basis(days,ts,flowPar,fdnames=fdnames)
  return(model)
}
# get models for observed and simulated time series
get_models <- function(obs.saved.file, sim.file, opt.zero, opt.fast, opt.new, loglam = seq(-4,2,0.25)){
  # read list of models for observed data
  obs.model.data <- readRDS(obs.saved.file)
  # read simulated time series
  data <- read.csv(sim.file, header = FALSE)
  # set date format
  data[,1] <- as.POSIXct(data[,1], format = '%Y-%m-%d')
  # remove additional day from leap years
  data <- data[format(data[,1], '%m-%d') != '02-29',]
  years.obs <- obs.model.data$best_model$fdnames$Years
  years.sim <- unique(format(data[,1], '%Y'))
  years     <- years.sim[years.sim %in% years.obs]
  # compute errors and coefficients
  obs.model   <- obs.model.data$best_model
  coefmat.sim <- numeric()
  eps0 <- rep(NA,length(years))
  eps1 <- rep(NA,length(years))
  lambda    <- obs.model.data$lambda
  flowbasis <- obs.model$basis
  if (opt.zero){
    opt.fast <- TRUE
  }
  if (opt.fast){
    if (length(lambda) > 1){
      opt.new <- TRUE
    }
    ts <- numeric()
    for (year in years){
      ts.annual  <- data[format(data[,1],'%Y') %in% year, ] # annual time-series
      ts <- cbind(ts,ts.annual[,2])
    }
    fdnames   <- list('Time','Years' = years, 'Variable')
    if (opt.zero){
      sim.model <- model_fit(ts,day.5,0,flowbasis,fdnames)$fd
    } else {
      if (opt.new){
        sim.model <- best_fit_fast(ts, loglam, 365, fdnames)$best_model$fd
      } else {
        sim.model <- model_fit(ts,day.5,lambda,flowbasis,fdnames)$fd
      }
    }
  } else {
    for (i in 1:length(years)){
      year         <- years[i]
      ts.annual    <- data[format(data[,1],'%Y') %in% year, ] # annual time-series
      lambda       <- obs.model.data$lambda[i]
      sim.model.i  <- model_fit(ts.annual[,2],day.5,lambda,flowbasis)
      coefmat.sim  <- cbind(coefmat.sim, coef(sim.model.i))
    }
    close(pb)
    fdnames   <- list('Time','Years' = years, 'Variable')
    sim.model <- fd(coefmat.sim, flowbasis, fdnames)
  }
  models <- list(sim.model = sim.model, obs.model = obs.model)
  return(models)
}
# get residuals for functions and first derivatives
get_res    <- function(models){
  sim.model  <- models$sim.model
  obs.model  <- models$obs.model
  deriv_sim.model <- deriv.fd(sim.model,1)
  deriv_obs.model <- deriv.fd(obs.model,1)
  diff_sq_fd  <- (sim.model - obs.model)^2
  Ddiff_sq_fd <- (deriv_sim.model - deriv_obs.model)^2
  const.basis <- create.constant.basis(c(0,365))
  epsilon0    <- inprod(const.basis,diff_sq_fd)
  epsilon1    <- inprod(const.basis,Ddiff_sq_fd)
  output <- list(eps0 = epsilon0, eps1 = epsilon1)
  return(output)
}
# plot results
plot_res   <- function(models, opt.plot, obs.saved.file){
  if (opt.plot){
    obs.model.data <- readRDS(obs.saved.file)
    days_list  <- obs.model.data$days_list
    sim.model  <- models$sim.model
    obs.model  <- models$obs.model
    deriv_sim.model <- deriv.fd(sim.model,1)
    deriv_obs.model <- deriv.fd(obs.model,1)
    years <- sim.model$fdnames$Years
    data_list <- obs.model.data$yobs
    for (i in 1:length(years)){
      year  <- years[i]
      if (is.null(days_list)){
        days <- day.5
        data <- data_list[,i]
      } else {
        days  <- days_list[[year]]
        data  <- data_list[[year]]
      }
      yobs  <- eval.fd(days,obs.model)[,i]
      ypred <- eval.fd(days,sim.model)[,i]
      Dyobs <- eval.fd(days,deriv_obs.model)[,i]
      Dysim <- eval.fd(days,deriv_sim.model)[,i]
      ylim <- range(yobs, ypred)
      plot(sim.model[i],ylim = c(0.8*ylim[1],1.1*ylim[2]), col = 'red', lty = 1, lwd = 2)
      par(new=TRUE)
      plot(obs.model[i],ylim = c(0.8*ylim[1],1.1*ylim[2]), lwd = 2)
      legend('topright', legend = c('Simulated','Observed'), col = c('red','black'), lty = c(1,1), lwd = c(2,2))
      par(new=FALSE)      
      ylim <- range(Dyobs, Dysim)
      plot(deriv_sim.model[i], ylim = 1.4*ylim, col ='red', lty = 1, lwd = 2)
      par(new=TRUE)
      plot(deriv_obs.model[i], ylim = 1.4*ylim)
      legend('topright', legend = c('Simulated','Observed'), col = c('red','black'), lty = c(1,1), lwd = c(2,2))
      par(new=FALSE)
    }
  }
}
# compute perfomance metric
metric     <- function(models, residuals, opt.plot){
  sim.model  <- models$sim.model
  obs.model  <- models$obs.model
  eps0 <- residuals$eps0
  eps1 <- residuals$eps1
  # t-test
  if (opt.plot){
    par()
    plot(mean(sim.model), col = 'red', lty = 1, lwd = 2)
    lines(mean(obs.model), lwd = 2)
    legend('topright', legend = c('Simulated','Observed'), col = c('red','black'), lty = c(1,1), lwd = c(2,2))
    par(mfrow = c(1,1))
  }
  t_out <- tperm.fd(obs.model,sim.model, plotres = opt.plot, argvals = day.5)
  # performance metric
  S <- mean(t_out$Tvals) # average
  of <- S*sqrt(sum(eps0+eps1))
  return(of)
}
# run routine to get performance metric 
get_metric <- function(obs.saved.file, sim.file, opt.plot, opt.zero, opt.fast, opt.new, loglam){
  # get the models
  models <- get_models(obs.saved.file, sim.file, opt.zero, opt.fast, opt.new, loglam)
  # compute residuals
  residuals <- get_res(models)
  # plot results
  plot_res(models, opt.plot, obs.saved.file)
  # obtain objective function
  of <- metric(models, residuals, opt.plot)
  return(of)
}