setwd("~/Google Drive/My Drive/cox_regression")
source("0_functions.R")
library(survival)
library(brms)
library(INLAjoint)
library(INLA)

#### Set parameters ####
n <- 200                # number of observations
p <- 2                  # number of covariates
beta <- c(0.5, -0.5)    # true coefficients for covariates
lambda <- 0.1           # scale parameter of Weibull (baseline hazard)
kappa <- 2              # shape parameter of Weibull (hazard increases with time)
censoring_rate <- 0.1   # desired rate of censoring

n_sim = 200
output_list = list()

#### log file code ####
# start new log file
write(paste0("Log file"," -- ",Sys.time()), file="log.txt")

for(i in 1:n_sim){
  # start time
  start_time <- Sys.time()
  #### Generate covariates ####
  X <- cbind( matrix(rnorm(n), nrow = n, ncol = 1), 
              matrix(runif(n,-1,1), nrow = n, ncol = 1)
              )
  colnames(X) <- paste0("X", 1:p)
  
  #### Generate survival times from Weibull baseline hazard ####
  # Generate sin(x) effect
  s = runif(n, 0, 2*pi)
  
  # Calculate linear predict
  eta <- X %*% beta + sin(s)
  # Using relationship between Weibull and Cox PH model: TT = (-log(U) / (lambda * exp(eta)))^(1 / kappa)
  U <- runif(n)
  TT <- (-log(U) / (lambda * exp(eta)))^(1 / kappa)
  
  #### Add censoring times ####
  C <- rexp(n, rate = censoring_rate)
  
  #### Observed times and event indicator ####
  time <- pmin(TT, C)
  status <- as.numeric(TT <= C)
  # print(mean(status))
  
  #### Create a data frame ####
  sim_data <- data.frame(time = time, status = status, X, s = s)
  
  # Inspect the data
  # head(sim_data)
  
  #### Fit Cox PH model ####
  cox_model <- coxph(Surv(time, status) ~ X1+X2+pspline(s, nterm=5), data = sim_data)
  GAM_fit = termplot(cox_model, term=3, se=TRUE, col.term=1, col.se=2, plot=F) # log hazard
  coxPH_summary = list(
    "coeff"=summary(cox_model),
    "GAM"=GAM_fit
    )
  
  #### Fit Cox-PG PH model ####
  coxPG1_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 n_print = 2500,burn_in = 1000,
                                 n_mcmc = 10000,epsilon = 1000,
                                 Ospline_variables = "s",
                                 thin = 10, calibration = F)
  coxPG1_summary = list(
    "coeff"=summary.coeff(coxPG1_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG1_model)$spline_summary,
    "time_grid"=coxPG1_model$rmb_splines$time_grid,
    "GAM"=summary.osplines(coxPG1_model,F)$spline_summary,
    "s_grid"=coxPG1_model$Osplines$s$domain + coxPG1_model$Ospline_means$s
  )
  
  coxPG2_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 n_print = 2500,burn_in = 1000,
                                 n_mcmc = 10000,
                                 Ospline_variables = "s",
                                 thin = 10, calibration = T)
  coxPG2_summary = list(
    "coeff"=summary.coeff(coxPG2_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG2_model)$spline_summary,
    "time_grid"=coxPG2_model$rmb_splines$time_grid,
    "GAM"=summary.osplines(coxPG2_model,F)$spline_summary,
    "s_grid"=coxPG2_model$Osplines$s$domain + coxPG2_model$Ospline_means$s
  )
  
  coxPG3_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 n_print = 2500,burn_in = 1000,
                                 n_mcmc = 10000, partitions = 10,
                                 Ospline_variables = "s",
                                 thin = 10, calibration = T)
  coxPG3_summary = list(
    "coeff"=summary.coeff(coxPG3_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG3_model)$spline_summary,
    "time_grid"=coxPG3_model$rmb_splines$time_grid,
    "GAM"=summary.osplines(coxPG3_model,F)$spline_summary,
    "s_grid"=coxPG3_model$Osplines$s$domain + coxPG3_model$Ospline_means$s
  )
  
  #### Fit INLA model ####
  INLA_model <- joint(formSurv = inla.surv(time = time, event = status) ~ X1+X2, 
                      basRisk = "rw1", dataSurv = sim_data)
  INLA_summary = summary(INLA_model)
  
  #### Fit HMC Bayesian Weibull model ####
  source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
  HMC_model <- brm(formula = bf(time | cens(status==0) ~ X1+X2+s(s,k=5)),
                   family = weibullPH,
                   stanvars = stanvars_weibullPH,
                   data = sim_data,
                   silent = 2,
                   refresh = 0,
                   warmup = 1000,
                   thin = 10,
                   chains = 1,       # Number of MCMC chains
                   iter = 11000,     # Number of iterations per chain
                   cores = 1         # Number of cores for parallel processing
  )
  time_plot = seq(min(time[status==1]),max(time[status==1]),length.out = 100)
  # Compute the baseline cumulative hazard function
  tmp = posterior_samples(HMC_model)
  tmp = exp(rep.col(tmp$b_Intercept,100)) * rep.row(time_plot,1000)^(rep.col(tmp$gamma,100))
  tmp = jointband_maps(log(tmp))
  tmp2 = conditional_smooths(HMC_model)
  HMC_summary = list(
    "coeff"=summary(HMC_model),
    "basline"=tmp,
    "time_grid"=time_plot,
    "GAM"=tmp2[[1]],
    "s_grid"=tmp2[[1]]$s
  )
  
  out = list(
    "coxPH"=coxPH_summary,
    "coxPG1"=coxPG1_summary,
    "coxPG2"=coxPG2_summary,
    "coxPG3"=coxPG3_summary,
    "INLA"=INLA_summary,
    "HMC"=HMC_summary
  )
  if(i%%1==0){
    cat("simulation #: ",i,"\n")
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    cat("time elapsed: ",time_taken,"\n")
    
    start_time <- Sys.time()
    
    #### log file code ####
    line=paste0("simulation #: ",i," -- ",Sys.time())
    write(line,file="log.txt",append=TRUE)
  }
  output_list[[i]] = out
}

saveRDS(output_list, file = "1_weibull_GAM.rds")