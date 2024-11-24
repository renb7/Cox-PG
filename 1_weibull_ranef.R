setwd("~/Google Drive/My Drive/cox_regression")
source("0_functions.R")
library(survival)
library(brms)
library(INLAjoint)
library(INLA)

#### Set parameters ####
n <- 200                # number of observations
n_clusters <- 25        # number of clusters
obs_per_cluster <- n / n_clusters  # observations per cluster
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
  # Assign each observation to a cluster
  cluster <- rep(1:n_clusters, each = obs_per_cluster)
  
  # Generate a random effect for each cluster
  b_cluster <- rnorm(n_clusters, mean = 0)
  # Assign each observation the random effect of its cluster
  b <- b_cluster[cluster]
  
  # Calculate linear predict
  eta <- X %*% beta + b
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
  sim_data <- data.frame(time = time, status = status, X, cluster = factor(cluster))
  cluster = model.matrix(~ as.factor(cluster) - 1)
  time_plot = seq(min(time[status==1]),max(time[status==1]),length.out = 100)
  
  # Inspect the data
  # head(sim_data)
  
  #### Fit Cox PH model ####
  cox_model <- coxph(Surv(time, status) ~ X1+X2+frailty.gaussian(cluster), data = sim_data)
  coxPH_summary = summary(cox_model)
  
  #### Fit Cox-PG PH model ####
  coxPG1_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 burn_in = 1000,
                                 n_mcmc = 11000,epsilon = 1000,
                                 Z_matrix_list = list("cluster"=cluster),
                                 thin = 10, calibration = F)
  coxPG1_summary = list(
    "coeff"=summary.coeff(coxPG1_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG1_model)$spline_summary,
    "time_grid"=coxPG1_model$rmb_splines$time_grid
  )
  
  coxPG2_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 burn_in = 1000,
                                 n_mcmc = 11000,
                                 Z_matrix_list = list("cluster"=cluster),
                                 thin = 10, calibration = T)
  coxPG2_summary = list(
    "coeff"=summary.coeff(coxPG2_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG2_model)$spline_summary,
    "time_grid"=coxPG2_model$rmb_splines$time_grid
  )
  
  coxPG3_model <- posterior_draw(Surv(time, status) ~ X1+X2, 
                                 data = sim_data,
                                 burn_in = 1000,
                                 n_mcmc = 11000, partitions = 10,
                                 Z_matrix_list = list("cluster"=cluster),
                                 thin = 10, calibration = T)
  coxPG3_summary = list(
    "coeff"=summary.coeff(coxPG3_model)$coeff_summary,
    "basline"=summary.baselineLogCumulativeHazard(coxPG3_model)$spline_summary,
    "time_grid"=coxPG3_model$rmb_splines$time_grid
  )
  
  #### Fit INLA model ####
  INLA_model <- joint(formSurv = inla.surv(time = time, event = status) ~ X1+X2+(1|id), 
                      basRisk = "rw1", id = "cluster", dataSurv = sim_data)
  # plot_out = plot(INLA_model)
  # plot_out = approx(plot_out$Baseline$data$ID,
  #                   exp(plot_out$Baseline$data$mean),
  #                   time_plot)$y*exp(-0.4164)
  # plot_out = cumulative_integral(time_plot,plot_out)
  INLA_summary = summary(INLA_model)
  
  #### Fit HMC Bayesian Weibull model ####
  source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
  HMC_model <- brm(formula = bf(time | cens(status==0) ~ X1+X2+(1|cluster)),
                   family =  weibullPH,
                   stanvars = stanvars_weibullPH,
                   data = sim_data, 
                   silent = 2,
                   refresh = 0,
                   warmup = 1000,
                   thin = 10,
                   chains = 1,       # Number of MCMC chains
                   iter = 11000      # Number of iterations per chain
                   )
  # Compute the baseline cumulative hazard function
  tmp = posterior_samples(HMC_model)
  tmp = exp(rep.col(tmp$b_Intercept,100)) * rep.row(time_plot,1000)^(rep.col(tmp$gamma,100))
  tmp = jointband_maps(log(tmp))
  HMC_summary = list(
    "coeff"=summary(HMC_model),
    "basline"=tmp,
    "time_grid"=time_plot
  )
  # output_list[[i]] = list(
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

saveRDS(output_list, file = "1_weibull_ranef.rds")