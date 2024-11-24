library(CVXR)
library(splines)
library(splines2)
library(BayesLogit)
library(survival)
library(relliptical)

cumulative_integral <- function(x, y) {
  # Ensure x and y are of equal length
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }
  
  # Initialize the cumulative integral array
  cum_integral <- numeric(length(x))
  cum_integral[1] <- 0  # Integral from the start to itself is zero
  
  # Compute the cumulative integral using the trapezoidal rule
  for (i in 2:length(x)) {
    dx <- x[i] - x[i - 1]
    cum_integral[i] <- cum_integral[i - 1] + (y[i] + y[i - 1]) * dx / 2
  }
  
  return(cum_integral)
}

symmetric_inv <- function(Q){
  sigma_new = solve(Q)
  sigma_new = (sigma_new+t(sigma_new))/2
  return(sigma_new)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

expit <- function(x){1/(1+exp(-x))}
logit <- function(p){log(p/(1-p))}

jointband_maps = function( MCMC_P,
                           alpha=c(0.05)){
  # Inputs
  #    'MCMC_P' - a B by K matrix containing MCMC samples. 
  #        B - number of MCMC iterations (MCMCspecs.B)
  #        K - number of function samples (number of columns in Y)
  #    'alpha' - a row vector containing significance levels at which 
  #        joint credible bands are to be returned. 
  
  # Outputs
  #   'MAPS' - multiplicity-adjusted probability scores. MAPS are
  #    truncated at either 0.5 or 0.001 for reducing computational burden.
  #   'upper_CI' - a length(alpha) by K matrix containing the upper bounds 
  #        of the joint credible bands. The first row corresponds to 
  #        the first level in alpha.       
  #   'lower_CI' - a length(alpha) by K matrix containing the lower bounds 
  #        of the joint credible bands. 
  B = dim(MCMC_P)[1]
  K = dim(MCMC_P)[2] 
  # sd_P = rep(NA,K)
  sd_P = apply(MCMC_P,MARGIN = 2,sd)
  mean_P = apply(MCMC_P,MARGIN = 2,mean)
  z_P = rep(NA,B)

  z_P.list = list()
  for(j in 1:B){
    # z_P is supremum over z-scores for each iteration
    # absolute value to make all z-scores positives, 
    # alpha does not get divided by 2
    z_P.list[[j]] = max( abs((MCMC_P[j,]-mean_P)/sd_P), na.rm = T )
  }
  
  
  z_P = unlist(z_P.list)
  
  levels = c(seq(0.5, 0.01, by = -0.01),
             seq(0.009, 0.0001, by = -0.001),
             seq(0.0009, 0.00001, by = -0.0001))
  cb1 = quantile(z_P,1-levels)
  MAPS = 0.5*rep(1,K)
  for(j in 2:length(levels)){
    temp_ind = (mean_P - cb1[j]*sd_P)*(mean_P + cb1[j]*sd_P)
    # (a-b)(a+b) = a^2 - b^2
    MAPS[temp_ind>0] = levels[j]
  }
  
  m = length(alpha)
  upper_CI = pracma::zeros(m,K)
  lower_CI = pracma::zeros(m,K)
  # use the empirical bounds based on z-scores
  # z_P is supremum over scores for each iteration
  cb2 = quantile(z_P,1-alpha);
  # i.e. quantile(z_P,1-0.05) should be about 1.96
  # absolute value to make all z-scores positives, 
  # alpha does not get divided by 2
  for(j in 1:m){
    upper_CI[j,] = mean_P + cb2[j]*sd_P;
    lower_CI[j,] = mean_P - cb2[j]*sd_P; 
  }
  return(list("MAPS"=MAPS,
              "upper_CI"=upper_CI,
              "lower_CI"=lower_CI)
  )
}

Osplines = function(x,numIntKnots=5,ngrid=1000,a=NULL,b=NULL){
  
  if(is.null(a) | is.null(b)){
    a = min(x)
    b = max(x)
  }
  
  xg <- seq(a,b,length=ngrid)
  
  intKnots <- quantile(unique(x),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
  
  names(intKnots) <- NULL
  B <- bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
  
  # dim(B) is number of intKnots + intercept + degree, 5+1+3=9
  # under default settings, should be 9
  
  # knots and boundaries
  # c(attr(B,"Boundary.knots")[1],attr(B,"knots"),attr(B,"Boundary.knots")[2])
  
  BTB <- crossprod(B,B) ; 
  # cat('intKnots=',intKnots,'\n')
  
  #### Derivative of B-spline basis
  formBprime <- function(x,a,b,intKnots)
  {
    allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
    K <- length(intKnots) 
    L <- 3*(K+8)
    Bd <- spline.des(allKnots,x,derivs=rep(1,length(x)),outer.ok=TRUE)$design     
    return(Bd)
  }
  
  formOmega <- function(a,b,intKnots)
  {
    allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
    K <- length(intKnots)
    L <- 3*(K+8)
    xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+rep(allKnots,each=3)[-c(1,2,L)])/2
    # dx in integral, used for numerical integration
    wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7) 
    # command to get spline design matrix
    Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design  
    Omega <- t(Bdd*wts)%*%Bdd     
    return(Omega)
  }
  
  # allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
  # B.des = spline.des(allKnots, x, outer.ok = T)$design
  # B.des = spline.des(allKnots, xg, outer.ok = T)$design
  
  Omega <- formOmega(a,b,intKnots)
  eigOmega <- eigen(Omega)
  
  # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:
  
  indsZ <- 1:(numIntKnots+2)
  UZ <- eigOmega$vectors[,indsZ] # Z_omega from paper
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ])) # diag term from paper
  
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- eigOmega$vectors[,indsX] # X_omega from paper
  
  # check to see that only two eigen values exist,
  # verify equation (17) from paper
  # cbind(B %*% UX,x) 
  # round(eigen(t(cbind(B %*% UX))%*%cbind(B %*% UX))$values,6)
  # round(eigen(t(cbind(B %*% UX,x))%*%cbind(B %*% UX,x))$values,6)
  # round(eigen(t(cbind(B %*% UX,x,1))%*%cbind(B %*% UX,x,1))$values,6)
  
  L <- cbind( UX, LZ )
  stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
  if (sum(stabCheck^2) > 1.0001*(numIntKnots+2)){
    print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")
  }
  
  X <- cbind(rep(1,length(x)),x)
  Z <- B%*%LZ # Demmler-Reinsch design matrix
  
  Bg <- bs(xg,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
  Zg <- Bg%*%LZ # Z matrix basis
  
  Bdg <- formBprime(xg,a,b,intKnots)
  Zdg <- Bdg%*%LZ # 1st derivative of Z matrix basis
  
  # Zres=rbind(Z,Zg,Zdg)
  # return(list(Zres=Zres,B=B,Omega=Omega,intKnots=intKnots))
  return(
    list(domain = xg,
         domain_obs = x,
         Z_all = Zg,
         DZ_all = Zdg,
         Z_obs = Z,
         B=B,Omega=Omega,intKnots=intKnots)
  )
}

rmb_splines <- function(surv.responses, partitions = 5, Z_g=NULL){
  all_times <- surv.responses[, "time"]
  max_time = max(all_times)
  all_times = all_times/(2*max_time)
  err = 0.0000001
  nt = 1000 # number of grid obs for splines vectors
  # Extract the status (event or censoring)
  status <- surv.responses[, "status"]
  if(!is.null(Z_g)){
    Du_obs_tmp = data.frame( "temp"=matrix(NA,nrow = length(all_times)) )
    u_obs_tmp = data.frame( "temp"=matrix(NA,nrow = length(all_times)) )
    Du_all_tmp = data.frame( "temp"=matrix(NA,nrow = nt) )
    u_all_tmp = data.frame( "temp"=matrix(NA,nrow = nt) )
    for(c in 1:ncol(Z_g)){
      # Get event times (status == 1) for each group
      event_times <- all_times[(status==1)*(Z_g[,c]==1)==1]
      
      J = min(length(unique(event_times)),partitions) # min on partitions, so each partition has a event
      time_seq = quantile(event_times, seq(0,1,length.out=J+1) )
      time_seq = sort(unique(time_seq)) # remove duplicates
      J = length(time_seq)-1
      time_seq[J+1] = time_seq[J+1] + err
      
      rmb_seq = seq( min(all_times), max(all_times), length.out=nt )
      rmb_seq = sort(unique(c(rmb_seq,time_seq))) # add partitions boundaries into continuous splines
      
      u_obs = data.frame( matrix(0,ncol=J,nrow = length(all_times)) )
      names(u_obs) = paste0("J",1:J)
      
      Du_obs = data.frame( matrix(0,ncol=J,nrow = length(all_times)) )
      names(Du_obs) = paste0("D",1:J)
      
      u_all = data.frame( matrix(0,ncol=J,nrow = length(rmb_seq)) )
      names(u_all) = paste0("J",1:J)
      
      Du_all = data.frame( matrix(0,ncol=J,nrow = length(rmb_seq)) )
      names(Du_all) = paste0("J",1:J)
      
      for(i in 2:(length(time_seq)) ){
        ind1 = rmb_seq < time_seq[i] & rmb_seq >= time_seq[i-1]
        Du_all[,i-1] = ind1*1
        
        ind2 = rmb_seq >= time_seq[i]
        u_all[ind1,i-1] = rmb_seq[ind1] - min(rmb_seq[ind1])
        u_all[ind2,i-1] = max(rmb_seq[ind1] - min(rmb_seq[ind1]))
        
        ind3 = all_times < time_seq[i] & all_times >= time_seq[i-1]
        Du_obs[,i-1] = ind3*1
        
        tmp = approx(x=rmb_seq, y=u_all[,i-1], xout = all_times)
        u_obs[,i-1] = tmp$y
      }
      # assign u_all on grid of nt
      rmb_seq = seq( min(all_times), max(all_times), length.out=nt )
      
      u_all = data.frame( matrix(0,ncol=J,nrow = nt) )
      names(u_all) = paste0("J",1:J)
      
      Du_all = data.frame( matrix(0,ncol=J,nrow = nt) )
      names(Du_all) = paste0("J",1:J)
      
      for(i in 2:(length(time_seq)) ){
        ind1 = rmb_seq < time_seq[i] & rmb_seq >= time_seq[i-1]
        Du_all[,i-1] = ind1*1
        
        ind2 = rmb_seq >= time_seq[i]
        u_all[ind1,i-1] = rmb_seq[ind1] - min(rmb_seq[ind1])
        u_all[ind2,i-1] = max(rmb_seq[ind1] - min(rmb_seq[ind1]))
      }
      
      # median centering
      u_obs = u_obs - rep.row(apply(u_all,2,median),nrow(u_obs))
      u_all = u_all - rep.row(apply(u_all,2,median),nrow(u_all))
      
      u_obs = u_obs*rep.col(Z_g[,c],ncol(u_obs))
      Du_obs = Du_obs*rep.col(Z_g[,c],ncol(Du_obs))
      
      colnames(Du_obs) = paste0( paste0("DZ_alpha_",1:J,"_",colnames(Z_g)[c]) )
      colnames(u_obs) = paste0( paste0("Z_alpha_",1:J,"_",colnames(Z_g)[c]) )
      
      colnames(Du_all) = paste0( paste0("DZ_alpha_",1:J,"_",colnames(Z_g)[c]) )
      colnames(u_all) = paste0( paste0("Z_alpha_",1:J,"_",colnames(Z_g)[c]) )
      
      # must use only deaths in density likelihood
      Du_obs = Du_obs*rep.col(status,ncol(Du_obs))
      
      Du_obs_tmp = cbind(Du_obs_tmp,Du_obs)
      u_obs_tmp = cbind(u_obs_tmp,u_obs)
      Du_all_tmp = cbind(Du_all_tmp,Du_all)
      u_all_tmp = cbind(u_all_tmp,u_all)
    }
    Du_obs_tmp$temp=NULL
    u_obs_tmp$temp=NULL
    Du_all_tmp$temp=NULL
    u_all_tmp$temp=NULL
    
    Du_obs = Du_obs_tmp
    u_obs = u_obs_tmp
    Du_all = Du_all_tmp
    u_all = u_all_tmp
  } else {
    # Get event times (status == 1)
    event_times <- all_times[status == 1]
    
    J = min(length(unique(event_times)),partitions) # min on partitions, so each partition has a event
    time_seq = quantile(event_times, seq(0,1,length.out=J+1) )
    time_seq = sort(unique(time_seq)) # remove duplicates
    J = length(time_seq)-1
    time_seq[J+1] = time_seq[J+1] + err
    
    rmb_seq = seq( min(all_times), max(all_times), length.out=nt )
    rmb_seq = sort(unique(c(rmb_seq,time_seq))) # add partitions boundaries into continuous splines
    
    u_obs = data.frame( matrix(0,ncol=J,nrow = length(all_times)) )
    names(u_obs) = paste0("J",1:J)
    
    Du_obs = data.frame( matrix(0,ncol=J,nrow = length(all_times)) )
    names(Du_obs) = paste0("D",1:J)
    
    u_all = data.frame( matrix(0,ncol=J,nrow = length(rmb_seq)) )
    names(u_all) = paste0("J",1:J)
    
    Du_all = data.frame( matrix(0,ncol=J,nrow = length(rmb_seq)) )
    names(Du_all) = paste0("J",1:J)
    
    for(i in 2:(length(time_seq)) ){
      ind1 = rmb_seq < time_seq[i] & rmb_seq >= time_seq[i-1]
      Du_all[,i-1] = ind1*1
      
      ind2 = rmb_seq >= time_seq[i]
      u_all[ind1,i-1] = rmb_seq[ind1] - min(rmb_seq[ind1])
      u_all[ind2,i-1] = max(rmb_seq[ind1] - min(rmb_seq[ind1]))
      
      ind3 = all_times < time_seq[i] & all_times >= time_seq[i-1]
      Du_obs[,i-1] = ind3*1
      
      tmp = approx(x=rmb_seq, y=u_all[,i-1], xout = all_times)
      u_obs[,i-1] = tmp$y
    }
    # assign u_all on grid of nt
    rmb_seq = seq( min(all_times), max(all_times), length.out=nt )
    
    u_all = data.frame( matrix(0,ncol=J,nrow = nt) )
    names(u_all) = paste0("J",1:J)
    
    Du_all = data.frame( matrix(0,ncol=J,nrow = nt) )
    names(Du_all) = paste0("J",1:J)
    
    for(i in 2:(length(time_seq)) ){
      ind1 = rmb_seq < time_seq[i] & rmb_seq >= time_seq[i-1]
      Du_all[,i-1] = ind1*1
      
      ind2 = rmb_seq >= time_seq[i]
      u_all[ind1,i-1] = rmb_seq[ind1] - min(rmb_seq[ind1])
      u_all[ind2,i-1] = max(rmb_seq[ind1] - min(rmb_seq[ind1]))
    }
    
    # median centering
    u_obs = u_obs - rep.row(apply(u_all,2,median),nrow(u_obs))
    u_all = u_all - rep.row(apply(u_all,2,median),nrow(u_all))
    
    colnames(Du_obs) = paste0( paste0("DZ_alpha_",1:J) )
    colnames(u_obs) = paste0( paste0("Z_alpha_",1:J) )
      
    colnames(Du_all) = paste0( paste0("DZ_alpha_",1:J) )
    colnames(u_all) = paste0( paste0("Z_alpha_",1:J) )
    
    # must use only deaths in density likelihood
    Du_obs = Du_obs*rep.col(status,ncol(Du_obs))
  }
  return(
    list(time_grid = 2*max_time*rmb_seq,
         time_obs = surv.responses[, "time"],
         u_all = u_all,
         Du_all = Du_all,
         u_obs = u_obs,
         Du_obs=Du_obs)
    )
}

M_matrix = function(formula, data, partitions = 10, Z_g=NULL){
  mf <- model.frame(formula=formula, data=data)
  X <- data.frame( model.matrix(formula, mf) )
  surv.responses <- model.response(mf)
  
  rmb_out = rmb_splines(surv.responses = surv.responses, partitions = partitions, Z_g=Z_g)
  
  # calculate intercept for each hazard group
  if(!is.null(Z_g)){
    intercept_new = matrix(X[,"X.Intercept."])
    for(c in 1:ncol(Z_g)){
      mat_tmp = data.frame( matrix(X[,"X.Intercept."]*rep.col(Z_g[,c],1)) )
      names(mat_tmp) = colnames( Z_g )[c]
      intercept_new = cbind(intercept_new,mat_tmp)
    }
    intercept_new = intercept_new[,-1]
    X = cbind(intercept_new,X[-1])
  }
  
  names(X)[names(X)=="X.Intercept."] = "intercept"
  
  M=cbind(rmb_out$u_obs,X)
  return(as.matrix(M))
}

initial_eta  <- function(y, Du_obs, M, w=NULL){
  # calculate constrained MLE
  if(is.null(w)){
    w=rep( 1, nrow(M) )
  }
  
  J = ncol(Du_obs)
  eta_cvx <- Variable(ncol(M))
  obj <- (-sum(y * w * M %*% eta_cvx) + sum( exp(M %*% eta_cvx)*w ) - sum( y*w*log(Du_obs%*%eta_cvx[c(1:J)]) ) )
  constraint1 <- eta_cvx[c(1:J)] >= 0
  
  prob <- Problem(Minimize(obj),constraints = list(constraint1))
  result <- CVXR::solve(prob, solver="SCS")
  eta_initial = result$getValue(eta_cvx)
  
  return(eta_initial)
}

log_MH_prob <- function(eta_new, eta_old, y, M, epsilon, w=NULL){
  if(is.null(w)){
    w=rep( 1, nrow(M) )
  }
  
  eta_ep = -log(epsilon)
  
  out1 = sum( exp(M %*% eta_old)*w )
  out1 = out1 + sum( -log( (1 + exp(M %*% eta_old + eta_ep))**(w*(y+epsilon)) ) )
  
  out2 = sum( exp(M %*% eta_new)*w )
  out2 = out2 + sum( -log( (1 + exp(M %*% eta_new + eta_ep))**(w*(y+epsilon)) ) )
  
  out = min(0, out1 - out2 )
  return( out )
}

posterior_draw <- function(formula, data, 
                           partitions = 5, epsilon=100,
                           eta_initial=NULL,
                           n_mcmc = 4000,
                           thin = 4,
                           hazard_group = NULL,
                           Ospline_variables = NULL,
                           Z_matrix_list = NULL,
                           weights = NULL,
                           n_print = 1000,
                           prior_sigma = NULL, prior_mean=NULL,
                           burn_in = 1000,
                           slice_burn_in = 100,
                           calibration = T,
                           verbose = F
                           ){
  # start time
  start_time <- Sys.time()
  
  eta_ep = -log(epsilon)
  mixed_model = !is.null(Z_matrix_list)
  # Center O-spline variables and create new formula
  Ospline_means = NULL
  if( !is.null(Ospline_variables) ){
    Ospline_means = list()
    for( i in Ospline_variables ){
      Ospline_means[[i]] = mean(data[,i])
      data[,i] = data[,i]-mean(data[,i])
    }
    formula = as.formula( paste0(deparse(formula),paste(paste0("+",Ospline_variables),collapse=" ")) )
    terms <- all.vars(formula[[3]])
    # Remove duplicate terms
    unique_terms <- unique(terms)
    # Rebuild the formula
    formula <- as.formula(paste(paste0(deparse(formula[[2]]),"~"), paste(unique_terms, collapse = " + ")))
  }
  # Get Ospline Z random effect matrices
  Z_Ospline_list = NULL
  if( !is.null(Ospline_variables) ){
    if( !mixed_model ){
      Z_matrix_list = list()
      mixed_model = T
    }
    Z_Ospline_list = list()
    for(z in 1:length(Ospline_variables)){
      O_out = Osplines(data[,Ospline_variables[z]])
      Z_matrix_list[[ Ospline_variables[z] ]] = O_out$Z_obs
      Z_Ospline_list[[ Ospline_variables[z] ]] = O_out
    }
  }
  # Get stratified hazard monotonic spline matrices
  hazard_group_factors = NULL
  if(!is.null(hazard_group)){
    data[,hazard_group] = as.factor(data[,hazard_group])
    intercept_ = data[,hazard_group]
    Z_g <- model.matrix(~ intercept_ - 1)
    colnames(Z_g) = levels(data[,hazard_group])
    hazard_group_factors = colnames(Z_g)
  } else {
    Z_g = NULL
  }
  # Rename Z random effect matrices columns and get dimensions
  MM_vec = c()
  Z_combined = c()
  col_names = c()
  if(mixed_model){
    for(z in 1:length(Z_matrix_list)){
      MM_vec = c( MM_vec, ncol(Z_matrix_list[[z]]) )
      Z_combined = cbind(Z_combined,Z_matrix_list[[z]])
      if( names(Z_matrix_list)[z] %in% Ospline_variables){
        col_names = c( col_names,paste0("Z_",names(Z_matrix_list)[z],"_", 1:ncol(Z_matrix_list[[z]]) ) )
      } else {
        col_names = c( col_names,paste0("Z_",colnames(Z_matrix_list[[z]])) ) 
      }
    }
    Z_combined = data.frame(Z_combined)
    colnames(Z_combined) = col_names
    MM = sum(MM_vec)
  } else {
    MM = 0
  }
  # Get design matrix
  M = M_matrix(formula, data, partitions = partitions, Z_g=Z_g)
  P = ncol(M) # get non random effect dimension
  M = cbind(M,Z_combined) # concat random effects
  
  if(is.null(prior_sigma)){
    A = diag(rep(1/1000000, P+MM))
  }else{
    A = diag(rep(1/1000000, P+MM))
    A[1:P,1:P] = symmetric_inv(prior_sigma)
  }
  
  if(is.null(prior_mean)){
    b_mu = rep(0,P+MM)
  }else{
    b_mu = rep(0,P+MM)
    b_mu[1+P] = prior_mean
  }
  # get death vector
  mf <- model.frame(formula=formula, data=data)
  surv.responses <- model.response(mf)
  y <- surv.responses[, "status"]*1
  # assign weights
  if(is.null(weights)){
    w=rep( 1, nrow(M) )
  } else {
    w=weights
  }
  # get monotonic splines
  rmb_out = rmb_splines(surv.responses = surv.responses, partitions = partitions, Z_g=Z_g)
  J = ncol(rmb_out$Du_obs)
  # get n_j beta parameters that accounts for weights
  nj = apply( rmb_out$Du_obs * rep.col(w,J), 2, sum)
  # get constrained MLE estimate
  eta0 = initial_eta(y, rmb_out$Du_obs, M, w=w)
  if(!is.null(eta_initial)){
    if(eta_initial=="zero"){
      eta0 = eta0*0
    } else {
      eta0 = matrix(eta_initial,ncol = 1)
    }
  }
  eta = eta0
  eta_new = t(eta)
  
  # set random effect precision tau
  tau_vec = c(rep(1/1000000, length(MM_vec)))
  tau_df = data.frame()
  M = as.matrix(M)
  eta_df = data.frame()
  accept_p = c()
  accept_rate = c()
  mcmc_num = c()
  for(iter in 1:(n_mcmc+burn_in)){
    
    if(iter%%n_print==0 && verbose){
      cat(iter,"\n")
      end_time <- Sys.time()
      time_taken <- end_time - start_time
      cat("time elapsed: ",time_taken,"\n")
    }
    
    # get bounds
    u_max = c()
    for(n_max in nj){
      v_tmp = rbeta(1,n_max,1)
      v_tmp = max(0,v_tmp)
      u_max = c(u_max,v_tmp)
    }
    v_max = u_max*eta[1:J]
    v_max = sapply(v_max, function(x) max(x,0) )
    u_plus = rep(1e4,J)
    # restrict lower and upper bound near MLE
    # v_max = rbind(v_max,eta0[1:J]*(1-0.5))
    # v_max = apply(v_max,2,max)
    # 
    # u_plus = rbind(u_plus,eta0[1:J]*(1+0.5))
    # u_plus = apply(u_plus,2,min)
    
    lower = as.matrix( c( v_max, rep(-Inf,(P+MM)-J) ) )
    upper = as.matrix( c( u_plus, rep(Inf,(P+MM)-J) ) )
    
    # get Gaussian moments
    omega = rpg(length(y), (y+epsilon)*w,
                (as.matrix(M)%*%eta) + eta_ep )
    omega = matrix(omega,nrow = length(y),ncol = 1)
    # for( inc in 1:length(y) ){
    #   omega[inc,] = rpg(1, (y[inc]+epsilon)*w[inc],
    #                     (as.matrix(M)%*%eta)[inc] + eta_ep )
    # }
    kappa = w*(y-epsilon)/2
    
    Q = t(M) %*% diag(omega[,1]) %*% M + A
    sigma_new = symmetric_inv(Q)
    mu = t(M) %*% ( kappa-omega*eta_ep ) + A %*% b_mu
    mu_new = sigma_new %*% mu
    
    # if(!matrixcalc::is.positive.definite(sigma_new)){
    #   sigma_new = sigma_new + diag(rep(1e-8,nrow(sigma_new)))
    # }
    
    # sample with bounds and Gaussian moments
    out = condMVNorm::condMVN(mean = mu_new,
                              sigma = sigma_new,
                              dependent.ind = c(1:J),
                              given.ind = c((J+1):(P+MM)),
                              X.given = c(eta[(J+1):(P+MM)]))
    mu_cond = t(t(out$condMean))
    sigma_cond = out$condVar
    sigma_cond = (t(sigma_cond) + sigma_cond)/2

    eta_new = rtelliptical(n = 1, mu = mu_cond[,1],
                           Sigma = sigma_cond,
                           lower = c(lower[1:J,1]),
                           upper = c(upper[1:J,1]),
                           dist = "Normal", burn.in = slice_burn_in,
                           thinning = 1)

    out = condMVNorm::condMVN(mean = mu_new,
                              sigma = sigma_new,
                              dependent.ind = c((J+1):(P+MM)),
                              given.ind = c(1:J),
                              X.given = c(eta_new))
    mu_cond = t(t(out$condMean))
    sigma_cond = out$condVar

    eta_cond = mu_cond + t(chol(sigma_cond)) %*% rnorm( ncol(sigma_cond) )
    eta_cond = t(eta_cond)

    eta_new = cbind(eta_new,eta_cond)
    
    # eta_new = rtelliptical(n = 1, mu = mu_new[,1],
    #                        Sigma = sigma_new, lower = c(lower[,1]), upper = c(upper[,1]),
    #                        dist = "Normal", nu = NULL, expr = NULL, gFun = NULL,
    #                        ginvFun = NULL, burn.in = slice_burn_in, thinning = 1)
    
    if(!calibration){ # debias towards Poisson process
      log_acc_p = 0
    } else {
      log_acc_p = log_MH_prob(t(eta_new), eta, y, M, epsilon, w = w)
    }
    accept_p = c(accept_p, exp(log_acc_p))
    acc = (log(runif(1))<log_acc_p)*1
    accept_rate = c(accept_rate, acc)
    
    if(acc){
      eta = eta_new[1,]
    }
    if(mixed_model){ # update precision
      A_tmp = rep(1/1000000, P+MM)
      MM_counter=0
      for(i in 1:length(MM_vec)){
        # print(MM_counter+P+(1:MM_vec[i]))
        eta_mixed = eta[MM_counter+P+(1:MM_vec[i])]
        tau_vec[i] = rgamma(1, 
                            shape = 1e-3 + MM_vec[i] / 2, 
                            rate = 1e-3 + sum(eta_mixed**2) / 2)
        A_tmp[MM_counter+P+(1:MM_vec[i])] = tau_vec[i]
        MM_counter = MM_counter + MM_vec[i]
      }
      A_tmp = diag(A_tmp)
      A_tmp[1:P,1:P] = A[1:P,1:P]
      A = A_tmp
    }
    if(iter%%thin==0 && iter > burn_in){
      eta_df = rbind(eta_df, t((eta)))
      tau_df = rbind(tau_df,tau_vec)
    }
  }
  # set column names
  colnames(eta_df) = colnames(M)
  colnames(tau_df) = names(Z_matrix_list)
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  if(verbose){
    cat("time taken: ",time_taken,"\n")
    cat("acceptance rate: ",mean(accept_rate),"\n")
  }
  return(
    list(
      posterior_samples = eta_df,
      tau_samples = tau_df,
      rmb_splines = rmb_out,
      y = y,
      M = M,
      initial_eta = eta0,
      formula = formula,
      accept_p = accept_p,
      accept_rate = accept_rate, 
      Ospline_variables = Ospline_variables,
      Osplines = Z_Ospline_list,
      Ospline_means = Ospline_means,
      hazard_group_factors = hazard_group_factors,
      Z_matrix_list = Z_matrix_list,
      P = P
    )
  )
}
summary.osplines <- function(fit, include_intercept=F){
  # get Z_osplines coefficient estimates
  coeff_output = list()
  coeff_summary = list()
  spline_output = list()
  spline_summary = list()
  for(z in fit$Ospline_variables){
    O_out = fit$Osplines[[z]]
    Z_matrix_names = c( z,paste0("Z_",z,"_", 1:ncol(fit$Z_matrix_list[[z]])) )
    # posterior draws
    coeff_draw = fit$posterior_samples[Z_matrix_names]
    # spline draw
    tmp = t( cbind( O_out$domain, O_out$Z_all) %*% t(coeff_draw) )
    # add in intercept
    if( is.null(fit$hazard_group_factors) ){
      intercept = rep.col(fit$posterior_samples[,"intercept"], ncol(tmp) )
      if(include_intercept){
        tmp = tmp + intercept
      }
    } else {
      intercept = t(t( rep(1/length(fit$hazard_group_factors),length(fit$hazard_group_factors)) ))
      intercept = as.matrix(fit$posterior_samples[,fit$hazard_group_factors]) %*% intercept
      intercept = rep.col(intercept, ncol(tmp) )
      if(include_intercept){
        tmp = tmp + intercept
      }
    }
    coeff_output[[z]] = coeff_draw
    spline_output[[z]] = tmp
    # get jointbands and credible intervals
    coeff_summary[[z]] = list(
      "mean" = apply(coeff_output[[z]],2,mean),
      "marginal_band" = apply(coeff_output[[z]],2, function(x) quantile(x,c(0.025,0.975))) 
    )
    spline_summary[[z]] = list( 
      "joint_band" = jointband_maps( spline_output[[z]] ),
      "mean" = apply(spline_output[[z]],2,mean),
      "marginal_band" = apply(spline_output[[z]],2, function(x) quantile(x,c(0.025,0.975))) 
      )
  }
  return(
    list(
      coeff_output = coeff_output,
      coeff_summary = coeff_summary,
      spline_output = spline_output,
      spline_summary = spline_summary
    )
  )
}
plot.osplines <- function(fit, term=1,...){
  out = summary.osplines(fit)
  i = names(out$spline_output)[term]
  jmap_out = out$spline_summary[[i]]$joint_band
  ylim = range(c(jmap_out$upper_CI,jmap_out$lower_CI))
  
  plot(fit$Osplines[[i]]$domain+fit$Ospline_means[[i]], 
       apply(out$spline_output[[i]],2,mean), ylim=ylim, type="l",
       ...)
  lines(fit$Osplines[[i]]$domain+fit$Ospline_means[[i]], 
        jmap_out$upper_CI, col="1", lty=2)
  lines(fit$Osplines[[i]]$domain+fit$Ospline_means[[i]], 
        jmap_out$lower_CI, col=1, lty=2)
}
summary.ranef <- function(fit){
  # get Z coefficient estimates
  coeff_output = list()
  coeff_summary = list()
  for(z in names(fit$Z_matrix_list)){
    if( !(z %in% fit$Ospline_variables) ){
      Z_matrix_names = paste0("Z_", colnames(fit$Z_matrix_list[[z]]) )
      coeff_output[[z]] = fit$posterior_samples[,Z_matrix_names]
      # get credible intervals
      coeff_summary[[z]] = list(
        "mean" = apply(coeff_output[[z]],2,mean),
        "marginal_band" = apply(coeff_output[[z]],2, function(x) quantile(x,c(0.025,0.975))) 
      )
    }
  }
  return(
    list(
      coeff_output = coeff_output,
      coeff_summary = coeff_summary
    )
  )
}
summary.baselineLogCumulativeHazard <- function(fit){
  # get monotonic coefficient estimates
  coeff_output = list()
  coeff_summary = list()
  spline_output = list()
  spline_summary = list()
  deriv_spline_output = list()
  if( !is.null(fit$hazard_group_factors) ){
    for(z in fit$hazard_group_factors){
      J = ncol(fit$rmb_splines$u_all)
      Z_matrix_names = paste0("Z_alpha_", 1:J, "_", z)
      Z_matrix_names = names(fit$rmb_splines$u_all)[names(fit$rmb_splines$u_all) %in% Z_matrix_names]
      J = length(Z_matrix_names)
      # posterior draws
      coeff_draw = fit$posterior_samples[Z_matrix_names]
      # spline draw
      tmp = t( as.matrix(fit$rmb_splines$u_all[Z_matrix_names]) %*% t(coeff_draw) )
      tmp2 = t( as.matrix(fit$rmb_splines$Du_all[paste0("D",Z_matrix_names)]) %*% t(coeff_draw) )
      # add in intercept
      intercept = rep.col(fit$posterior_samples[,z], ncol(tmp) )
      tmp = tmp + intercept
      coeff_output[[z]] = fit$posterior_samples[c(Z_matrix_names,z)]
      spline_output[[z]] = tmp
      deriv_spline_output[[z]] = tmp2
      # get jointbands and credible intervals
      coeff_summary[[z]] = list(
        "mean" = apply(coeff_output[[z]],2,mean),
        "marginal_band" = apply(coeff_output[[z]],2, function(x) quantile(x,c(0.025,0.975))) 
      )
      spline_summary[[z]] = list( 
        "joint_band" = jointband_maps( spline_output[[z]] ),
        "mean" = apply(spline_output[[z]],2,mean),
        "marginal_band" = apply(spline_output[[z]],2, function(x) quantile(x,c(0.025,0.975))) 
      )
    }
  } else {
    J = ncol(fit$rmb_splines$u_all)
    Z_matrix_names = paste0("Z_alpha_", 1:J)
    # posterior draws
    coeff_draw = fit$posterior_samples[Z_matrix_names]
    # spline draw
    tmp = t( as.matrix(fit$rmb_splines$u_all) %*% t(coeff_draw) )
    tmp2 = t( as.matrix(fit$rmb_splines$Du_all) %*% t(coeff_draw) )
    # add in intercept
    intercept = rep.col(fit$posterior_samples[,"intercept"], ncol(tmp) )
    # print(dim(intercept))
    tmp = tmp + intercept
    coeff_output[["base"]] = fit$posterior_samples[c(Z_matrix_names,"intercept")]
    spline_output[["base"]] = tmp
    deriv_spline_output[["base"]] = tmp2
    # get jointbands and credible intervals
    coeff_summary[["base"]] = list(
      "mean" = apply(coeff_output[["base"]],2,mean),
      "marginal_band" = apply(coeff_output[["base"]],2, function(x) quantile(x,c(0.025,0.975))) 
    )
    spline_summary[["base"]] = list( 
      "joint_band" = jointband_maps( spline_output[["base"]] ),
      "mean" = apply(spline_output[["base"]],2,mean),
      "marginal_band" = apply(spline_output[["base"]],2, function(x) quantile(x,c(0.025,0.975))) 
    )
  }
  return(
    list(
      coeff_output = coeff_output,
      coeff_summary = coeff_summary,
      spline_output = spline_output,
      spline_summary = spline_summary,
      deriv_spline_output = deriv_spline_output
    )
  )
}
summary.coeff <- function(fit){
  # get coefficient estimates
  coeff_output = list()
  coeff_summary = list()
  
  J = ncol(fit$rmb_splines$u_obs)
  a = max( length( fit$hazard_group_factors ), 1 ) # index adjustment
  P = fit$P
  # posterior draws
  coeff_draw = fit$posterior_samples[(J+a+1):P]
  coeff_summary = list(
    "mean" = apply(coeff_draw,2,mean),
    "marginal_band" = apply(coeff_draw,2, function(x) quantile(x,c(0.025,0.975))) 
  )
  
  return(
    list(
      coeff_output = coeff_draw,
      coeff_summary = coeff_summary
    )
  )
}
summary.coeffCombination <- function(fit, coeff_list){
  # get coefficient estimates
  spline_eff = rep(0,nrow(fit$posterior_samples))
  out_new = summary.osplines(fit, include_intercept = F)
  for(z in names(out_new$spline_output)){
    Zdomain = fit$Osplines[[z]]$domain + fit$Ospline_means[[z]]
    out_tmp = out_new$spline_output[[z]]
    out_tmp = out_tmp[,which.min( abs(Zdomain-coeff_list[[z]]) )]
    spline_eff = spline_eff + out_tmp
  }
  
  
  J = ncol(fit$rmb_splines$u_obs)
  a = max( length( fit$hazard_group_factors ), 1 ) # index adjustment
  P = fit$P
  # posterior draws
  coeff_draw = fit$posterior_samples[(J+a+1):P]
  for(i in names(coeff_draw)){
    if(i %in% names(coeff_list)){
      if( !(i %in% names(out_new$spline_output)) ){
        spline_eff = spline_eff + coeff_draw[,i]*coeff_list[[i]]
      }
    }
  }
  
  return(
    list(
      eff_output = spline_eff
    )
  )
}
plot.baselineLogCumulativeHazard <- function(fit){
  out = summary.baselineLogCumulativeHazard(fit)
  ylim = c()
  for( i in names(out$coeff_output) ){
    jmap_out = out$spline_summary[[i]]$joint_band
    ylim = c(ylim, range(c(jmap_out$upper_CI,jmap_out$lower_CI)) )
  }
  ylim = range(ylim)
  
  col=1
  plot(fit$rmb_splines$time_grid, 
       apply(out$spline_output[[1]],2,mean), 
       col=col, ylim=ylim, type="l")
  for( i in names(out$coeff_output) ){
    jmap_out = out$spline_summary[[i]]$joint_band
    lines(fit$rmb_splines$time_grid, 
          apply(out$spline_output[[i]],2,mean), col=col)
    lines(fit$rmb_splines$time_grid, 
          jmap_out$upper_CI, col=col, lty=2)
    lines(fit$rmb_splines$time_grid, 
          jmap_out$lower_CI, col=col, lty=2)
    col=col+1
  }
}
plot.baselineSurvival <- function(fit, coeff_list = NULL){
  out = summary.baselineLogCumulativeHazard(fit)
  ylim = c()
  jmap_list = list()
  mean_list = list()
  if(!is.null(coeff_list)){
    out2 = summary.coeffCombination(fit, coeff_list)
    for( i in names(out$coeff_output) ){
      jmap_out = out$spline_output[[i]] + rep.col( out2$eff_output, ncol(out$spline_output[[i]]) )
      mean_list[[i]] = apply(jmap_out,2,mean)
      jmap_out = jointband_maps( jmap_out )
      jmap_list[[i]] = jmap_out
      ylim = c(ylim, range(c(jmap_out$upper_CI,jmap_out$lower_CI)) )
    }
    ylim = range(exp(-exp(ylim)))
  } else {
    for( i in names(out$coeff_output) ){
      jmap_out = out$spline_output[[i]]
      mean_list[[i]] = apply(jmap_out,2,mean)
      jmap_out = jointband_maps( jmap_out )
      jmap_list[[i]] = jmap_out
      ylim = c(ylim, range(c(jmap_out$upper_CI,jmap_out$lower_CI)) )
    }
    ylim = range(exp(-exp(ylim)))
  }
  
  col=1
  plot(fit$rmb_splines$time_grid, 
       exp(-exp(mean_list[[1]])), 
       col=col, ylim=ylim, type="l")
  for( i in names(mean_list) ){
    jmap_out = jmap_list[[i]]
    lines(fit$rmb_splines$time_grid, 
          exp(-exp(mean_list[[i]])), col=col)
    lines(fit$rmb_splines$time_grid, 
          exp(-exp(jmap_out$upper_CI)), col=col, lty=2)
    lines(fit$rmb_splines$time_grid, 
          exp(-exp(jmap_out$lower_CI)), col=col, lty=2)
    col=col+1
  }
}
# plot.baselineHazard <- function(fit){
#   out = summary.baselineLogCumulativeHazard(fit)
#   out_new = list()
#   jmap_new = list()
#   ylim = c()
#   for( i in names(out$coeff_output) ){
#     out_new[[i]] = exp(out$spline_output[[i]]) * out$deriv_spline_output[[i]]
#     jmap_new[[i]] = jointband_maps( out_new[[i]] )
#     ylim = c(ylim, range(c(jmap_new[[i]]$upper_CI,jmap_new[[i]]$lower_CI)) )
#   }
#   ylim = range(ylim)
# 
#   col=1
#   plot(fit$rmb_splines$time_grid,
#        apply(out_new[[1]],2,mean),
#        col=col, ylim=ylim, type="l")
#   for( i in names(out$coeff_output) ){
#     jmap_out = jmap_new[[i]]
#     lines(fit$rmb_splines$time_grid,
#           apply(out_new[[i]],2,mean), col=col)
#     lines(fit$rmb_splines$time_grid,
#           jmap_out$upper_CI, col=col, lty=2)
#     lines(fit$rmb_splines$time_grid,
#           jmap_out$lower_CI, col=col, lty=2)
#     col=col+1
#   }
# }
