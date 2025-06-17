# replication_functions.R
#
# -----------------------------------------------------------------------------
# Load required packages
# -----------------------------------------------------------------------------


library(here)           # for project-relative paths
library(dirmult)        # rdirichlet for Dirichlet sampling
library(pracma)         # randortho and pinv
library(MASS)           # mvrnorm
library(abind)          # abind for array binding
library(foreach)        # foreach loops
library(doParallel)     # parallel backend for foreach
library(parallel)       # detectCores, makeCluster
library(stats)          # rnorm, prcomp, etc.
library(haven)          # read_dta
library(dplyr)


# Then, the potential outcomes are defined by Y(0) = Z^T \beta + \epsilon
# where Z are the principal components of X and \epsilon ~ N(0, 1). We use
# an additive treatment effect such that Y(1) = Y(0) + \tau
generate_data <- function(n, d, tau, conc, beta.choice, p=0.5) {
  
  P = randortho(n=d) #Random orthogonal matrix
  eigs = sort(rdirichlet(n=1, alpha=rep(conc, d)) * d, decreasing=TRUE)
  mean.vector = rep(0, d) #Mean vector
  cov.matrix = P %*% diag(as.list(eigs)) %*% t(P)
  
  x = mvrnorm(n=n, mean.vector, cov.matrix)
  
  # For the case that beta is a function of the eigenvalues
  if(beta.choice == 'uniform') {
    beta.z = as.matrix(rep(1, d))
  } else if(beta.choice == 'sqrt') {
    beta.z = as.matrix(sort(sqrt(eigs), decreasing=TRUE))
  } else {
    beta.z = as.matrix(sort(sqrt(eigs), decreasing=FALSE))
  }
  
  svd = prcomp(x, center=TRUE, scale=TRUE)
  pcs = svd$x
  y0 = rnorm(n, pcs%*%beta.z, 1)
  y1 = y0 + tau
  
  data = data.frame(cbind(scale(data.frame(x),
                                center=TRUE, scale=TRUE),
                          y0=y0, y1=y1))
  
  pca.data = data.frame(cbind(pcs, y0=y0, y1=y1))
  
  return(list(data, pca.data))
}


# Finds the number of principal components to include as decided by the
# Kaiser rule. We include all components such that their variance explained
# is greater than the average amount of variance explained.
get_kaiser_dim = function(eigs, d) {
  
  if (d == 1) {
    return(d)
  }
  
  var_explained = eigs / sum(eigs)
  avg_explained = mean(var_explained)
  
  return(min(max(which(var_explained > avg_explained)), d))
}


# Returns the vector (v_a(k), v_a(k), ..., v_a(k), 1, ..., 1)
pca_nu = function(pca_k, d, alpha=0.05) {
  
  pca_a = qchisq(alpha, df=pca_k)
  pca_nu = pchisq(pca_a, df=(pca_k+2)) / pchisq(pca_a, df=pca_k)
  
  return(c(rep(pca_nu, pca_k), rep(1, d - pca_k)))
  
}


# Creates randomized treatment vector and returns dataset
randomize = function(xyData, p=0.5){
  n = nrow(xyData)
  
  #half to treatment, half to control
  treatment = c(rep(1, n*p), rep(0, n*(1-p)))
  treatment = sample(treatment)
  
  xData = subset(xyData, select = -c(y0, y1))
  #if assigned control, get y0; otherwise get y1
  y = ifelse(treatment == 0, xyData$y0, xyData$y1)
  
  return(cbind(xData, y, treatment))
}


# Difference in means calculation
getMeanDiff = function(data){
  
  treatmentData = subset(data, treatment == 1)
  controlData = subset(data, treatment == 0)
  
  return(mean(treatmentData$y) - mean(controlData$y))
}


# Special case of get_weighted_k_simulated_eigs where we are using v_a(k)
# instead of \nu_{j, \eta}(k). Because v_a(k) has a closed-form expression
# this function is significantly faster and does not rely on Monte-Carlo
# simulation.
get_weighted_k_mahalanobis = function(xyData, pca.data, eigs, d, alpha=0.05) {
  
  # If d = 1 then no selection is necessary
  if(d == 1) {
    return(d)
  }
  
  full_a = qchisq(alpha, df=d)
  full_nu = pchisq(full_a, df=(d+2)) / pchisq(full_a, df=d)
  
  weighted_k = d # Initialize choice of k
  
  lb = 0 # Setting zero as an initial lower-bound
  for(j in 1:(d-1)) {
    
    # The following code obtains the variance reduction factors
    # v_a(d) and v_a(k) for the full and reduced Mahalanobis-based models,
    # respectively
    reduced_a = qchisq(alpha, df=j)
    reduced_nu = pchisq(reduced_a, df=(j+2)) / pchisq(reduced_a, df=j)
    
    
    # Cost-benefit calculation weighted by the eigenvalues of \Sigma
    benefit = sum(eigs[1:j] * (full_nu - reduced_nu))
    cost = sum(eigs[(j+1):d] * (1 - full_nu))
    
    if (benefit - cost >= lb) {
      lb = benefit - cost
      weighted_k = j
    }
  }
  
  return(weighted_k)
}

# Used to obtain a valid randomization according to the threshold a.
# We return z which will be used to calculate the variance reduction factor \nu
get_acceptable_rand <- function(eigs, a, d) {
  ub = as.numeric('inf')
  while (ub > a) {
    z = rchisq(n=d, df=1)
    ub = sum(eigs * z)
  }
  return(z)
}


# Averages the value of Z^2_j | \sum^d_{j=1} \eta_j Z^2_j
calculate_nu <- function(eigs, a, d, sims=1000) {
  z_list = lapply(1:sims, function(x) get_acceptable_rand(eigs, a, d))
  return(colMeans(do.call(rbind, z_list)))
}


# Used to calculate the threshold a for a given set of eigenvalues. We
# simulate \sum^d_{j=1} \eta_j Z^2_j sims many times, then take the
# desired quantile of the distribution
calculate_a <- function(eigs, d, alpha=0.05, sims=1000) {
  qf_dist = sapply(1:sims, function(x) sum(eigs * rchisq(n=d, df=1)))
  return(quantile(qf_dist, probs=alpha))
  
}


# Monte-Carlo approximation of \nu_{j, \eta}
simulate_nu = function(eigs, d, alpha=0.05, a_sims=10000, nu_sims=1000) {
  simulated_a = calculate_a(eigs=eigs, d=d, alpha=alpha, sims=a_sims)
  return(calculate_nu(eigs=eigs, a=simulated_a, d=d, sims=nu_sims))
}


# Function for simulating Figure 1
simulate_nu_df = function(n, d, tau, conc, beta, alpha=0.05, p=0.5) {
  
  # Note that the choice of beta.z is not important for this simulation
  sim.data = generate_data(n=n, d=d, tau=tau, conc=conc, beta='uniform')
  
  xyData = sim.data[[1]] # Simulated dataframe
  pca.df = sim.data[[2]] 
  
  covariates = as.matrix(xyData[, 1:d]) #Isolate X_1, ..., X_d
  eigs = n*p*(1 - p)*eig(cov(covariates)) # Eigenvalues of Cov(X)
  
  kaiser_k = get_kaiser_dim(eigs=eigs, d=d)
  kaiser_nu = pca_nu(pca_k=kaiser_k, d=d, alpha=alpha) #PCA Rerandomization using Kaiser Rule
  
  weighted_k = get_weighted_k_mahalanobis(xyData=xyData, pca.data=pca.df, eigs=eigs, d=d, alpha=alpha)
  weighted_nu = pca_nu(pca_k=weighted_k, d=d, alpha=alpha) #PCA Rerandomization using Weighted Eigenvalue Rule
  
  # Value of \nu under Mahalanobis Rerandomization
  mdr_a = qchisq(alpha, df=d)
  mdr_nu_single = pchisq(mdr_a, df=(d+2)) / pchisq(mdr_a, df=d)
  mdr_nu = rep(mdr_nu_single, d)
  
  edr_nu = simulate_nu(eigs=eigs, d=d, alpha=alpha, a_sims=10000, nu_sims=1000) #Euclidean Rerandomization
  sq_edr_nu = simulate_nu(eigs=(eigs ** 2), d=d, alpha=alpha, a_sims=10000, nu_sims=1000) #Squared Euclidean Rerandomization
  
  nu.df <- data.frame(Mahalanobis=mdr_nu, Euclidean=edr_nu,
                      Squared_Euclidean=sq_edr_nu,
                      Kaiser=kaiser_nu,
                      Weighted=weighted_nu, idx=1:d)
  return(nu.df)
}


# Used for Monte-Carlo simulation of the threshold a
get_threshold <- function(data, inv_cov, metric, sims=1000, alpha=0.05, p=0.5, ...) {
  
  return(quantile(sapply(1:sims,
                         function(i) metric(data=randomize(data, p),
                                            inv_cov=inv_cov, ...)),
                  probs=alpha)) 
}


# Used for Monte-Carlo simulation of the threshold a
get_threshold_par <- function(data, inv_cov, metric, sims=1000, alpha=0.05, p=0.5, ...) {
  
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  vals <- foreach(i=1:sims, .export = c("randomize", "getCovs")) %dopar% {
    metric(data=randomize(data, p), inv_cov=inv_cov, ...)
  }
  
  
  stopCluster(cl)
  
  return(quantile(unlist(vals), probs=alpha)) 
}


# Obtains just the covariates X_1, ..., X_d from a randomized data set
getCovs = function(data){
  return(subset(data, select = -c(treatment, y)))
}


# Mahalanobis Distance
get_MD <- function(data, inv_cov, p=0.5) {
  
  n = nrow(data)
  treatmentData = getCovs(subset(data, treatment == 1))
  controlData = getCovs(subset(data, treatment == 0))
  
  diff = sqrt(n) * (colMeans(treatmentData) - colMeans(controlData))
  
  return(t(diff) %*% inv_cov %*% diff)
}


# "Squared" Euclidean Distance
get_SED <- function(data, inv_cov, p=0.5) {
  
  n = nrow(data)
  treatmentData = getCovs(subset(data, treatment == 1))
  controlData = getCovs(subset(data, treatment == 0))
  
  n = nrow(data)
  covX = cov(getCovs(data)) / (p * (1 - p))
  
  diff = sqrt(n) * (colMeans(treatmentData) - colMeans(controlData))
  
  return(t(diff) %*% covX %*% diff)
}


# Euclidean Distance
get_ED <- function(data, inv_cov, p=0.5) {
  
  n = nrow(data)
  treatmentData = getCovs(subset(data, treatment == 1))
  controlData = getCovs(subset(data, treatment == 0))
  
  diff = sqrt(n) * (colMeans(treatmentData) - colMeans(controlData))
  
  return(t(diff) %*% diff)
}


# Oracle distance
get_opt_distance <- function(data, inv_cov, beta, p=0.5) {
  
  n = nrow(data)
  treatmentData = getCovs(subset(data, treatment == 1))
  controlData = getCovs(subset(data, treatment == 0))
  
  diff = sqrt(n) * (colMeans(treatmentData) - colMeans(controlData))
  
  return(t(diff) %*% (beta %*% t(beta)) %*% diff)
}


# Obtains an acceptable randomization where Q_A(x) <= a
get_rerand <- function(data, inv_cov, threshold, metric, p, ...) {
  
  ub = as.numeric('inf')
  while (ub > threshold) {
    randomization = randomize(data, p=p)
    ub = metric(randomization, inv_cov, ...)
    
  }
  return(randomization)
}


simulate_tauhat <- function(n, d, tau, conc, beta.choice, alpha, p=0.5, a_sims=100) {
  sim.data = generate_data(n=n, d=d, tau=tau, conc=conc, beta.choice=beta.choice)
  
  xyData = sim.data[[1]] # Simulated dataframe
  pca.df = sim.data[[2]] 
  
  n = nrow(xyData)
  
  covariates = as.matrix(xyData[, 1:d]) #Isolate X_1, ..., X_d
  eigs = n*p*(1 - p)*eig(cov(covariates)) # Eigenvalues of Cov(X)
  inv_cov_matrix = pinv(cov(covariates)/(p*(1-p)))
  
  s.y0.x = t(covariates) %*% as.matrix(xyData$y0 - mean(xyData$y0)) / (n-1)
  s.y1.x = t(covariates) %*% as.matrix(xyData$y1 - mean(xyData$y1)) / (n-1)
  v.x.tau = (1 / p) * s.y1.x + (1 / (1-p)) * s.y0.x
  
  beta = inv_cov_matrix %*% v.x.tau
  
  kaiser_k = get_kaiser_dim(eigs=eigs, d=d)
  kaiser_nu = pca_nu(pca_k=kaiser_k, d=d, alpha=alpha) #PCA Rerandomization using Kaiser Rule
  kaiser_df = data.frame(cbind(pca.df[,1:kaiser_k], y0=xyData$y0, y1=xyData$y1))
  kaiser_inv = solve(cov(as.matrix(kaiser_df[, 1:kaiser_k])))
  
  weighted_k = get_weighted_k_mahalanobis(xyData=xyData, pca.data=pca.df, eigs=eigs, d=d, alpha=alpha)
  weighted_nu = pca_nu(pca_k=weighted_k, d=d, alpha=alpha) #PCA Rerandomization using Weighted Eigenvalue Rule
  weighted_df = data.frame(cbind(pca.df[,1:weighted_k], y0=xyData$y0, y1=xyData$y1))
  weighted_inv = pinv(cov(as.matrix(weighted_df[, 1:weighted_k])))
  
  md_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_MD, sims=a_sims, alpha=alpha, p=p)
  ed_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_ED, sims=a_sims, alpha=alpha, p=p)
  sed_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_SED, sims=a_sims, alpha=alpha, p=p)
  kaiser_a <- get_threshold(data=kaiser_df, inv_cov=kaiser_inv, metric=get_MD, sims=a_sims, alpha=alpha, p=p)
  weighted_a <- get_threshold(data=weighted_df, inv_cov=weighted_inv, metric=get_MD, sims=a_sims, alpha=alpha, p=p)
  opt_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_opt_distance,
                         sims=a_sims, alpha=alpha, p=p, beta=beta)
  
  
  complete_var_rerands <- randomize(xyData, p=p)
  md_rerands <- get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=md_a, metric=get_MD, p=p)
  ed_rerands <- get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=ed_a, metric=get_ED, p=p)
  sed_rerands <- get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=sed_a, metric=get_SED, p=p)
  kaiser_rerands <- get_rerand(data=kaiser_df, inv_cov=kaiser_inv, threshold=kaiser_a, metric=get_MD, p=p)
  weighted_rerands <- get_rerand(data=weighted_df, inv_cov=weighted_inv, threshold=weighted_a, metric=get_MD, p=p)
  opt_rerands <- get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=opt_a, metric=get_opt_distance, p=p, beta=beta)
  
  complete.var <- getMeanDiff(complete_var_rerands)
  md.var = getMeanDiff(md_rerands)
  ed.var = getMeanDiff(ed_rerands)
  sed.var = getMeanDiff(sed_rerands)
  kaiser.var = getMeanDiff(kaiser_rerands)
  weighted.var = getMeanDiff(weighted_rerands)
  opt.var = getMeanDiff(opt_rerands)
  
  
  var.df = data.frame(cbind('Complete Randomization'=complete.var,
                            'Mahalanobis'=md.var,
                            'Euclidean'=ed.var,
                            'Sq_Euclidean'=sed.var,
                            'Kaiser'=kaiser.var,
                            'Weighted'=weighted.var,
                            'Optimal'=opt.var))
  
  return(var.df)
}  


simulate_var_study <- function(n, d, tau, conc, p, alpha, rerands, a_sims) {
  
  
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  vals.uniform <- foreach(i=1:rerands, .combine = rbind,
                          .export = c("simulate_tauhat", "get_rerand",
                                      "getCovs", "getMeanDiff", "randomize",
                                      "get_MD", "get_ED", "get_SED",
                                      "get_opt_distance", "generate_data",
                                      "get_threshold", "get_kaiser_dim",
                                      "pca_nu", "get_weighted_k_mahalanobis"),
                          .packages=c('MASS', 'pracma', 'dirmult',
                                      'magrittr', 'dplyr')) %dopar% {
                                        simulate_tauhat(n=n, d=d, tau=tau,
                                                        conc=conc,
                                                        beta.choice='uniform',
                                                        alpha=alpha, p=p,
                                                        a_sims=a_sims)}
  
  
  
  vals.sqrt <- foreach(i=1:rerands, .combine = rbind,
                       .export = c("simulate_tauhat", "get_rerand", "getCovs",
                                   "getMeanDiff", "randomize", "get_MD",
                                   "get_ED", "get_SED", "get_opt_distance",
                                   "generate_data","get_threshold",
                                   "get_kaiser_dim", "pca_nu",
                                   "get_weighted_k_mahalanobis"),
                       .packages=c('MASS', 'pracma', 'dirmult',
                                   'magrittr', 'dplyr')) %dopar% {
                                     simulate_tauhat(n=n, d=d, tau=tau,
                                                     conc=conc,
                                                     beta.choice='sqrt',
                                                     alpha=alpha, p=p,
                                                     a_sims=a_sims)}
  
  
  
  vals.inv <- foreach(i=1:rerands, .combine = rbind,
                      .export = c("simulate_tauhat","get_rerand",
                                  "getCovs", "getMeanDiff", "randomize",
                                  "get_MD", "get_ED", "get_SED",
                                  "get_opt_distance", "generate_data",
                                  "get_threshold", "get_kaiser_dim", "pca_nu",
                                  "get_weighted_k_mahalanobis"),
                      .packages=c('MASS', 'pracma', 'dirmult',
                                  'magrittr', 'dplyr')) %dopar% {
                                    simulate_tauhat(n=n, d=d, tau=tau,
                                                    conc=conc,
                                                    beta.choice='inv',
                                                    alpha=alpha, p=p,
                                                    a_sims=a_sims)}
  
  stopCluster(cl)
  
  
  uniform.sd = data.frame(as.list(apply(vals.uniform, MARGIN=2, sd)))
  var.uniform <- cbind(uniform.sd, beta='uniform', d=d, Concentration=conc)
  
  sqrt.sd <- data.frame(as.list(apply(vals.sqrt, MARGIN=2, sd)))
  var.sqrt <- cbind(sqrt.sd, beta='sqrt', d=d, Concentration=conc)
  
  inv.sd <- data.frame(as.list(apply(vals.inv, MARGIN=2, sd)))
  var.inv <- cbind(inv.sd, beta='inv', d=d, Concentration=conc)
  
  out_file <-here("output", "calculations", sprintf("var_%s_%s_alpha_01_a2500.csv", d, conc))
  
  write.csv(rbind(var.sqrt, var.uniform, var.inv), out_file, row.names = FALSE)
  
  
  return(rbind(var.sqrt, var.uniform, var.inv))
  
}


get_weighted_k_sim = function(eigs, d, alpha=0.05, a_sims=10000, nu_sims=10000) {
  
  # If d = 1 then no selection is necessary
  if(d == 1) {
    return(1)
  }
  
  weighted_k = d
  full_nu = simulate_nu(eigs=eigs, d=d, alpha=alpha, a_sims=a_sims, nu_sims=nu_sims)
  lb = 0 # Setting zero as an initial lower-bound
  for(j in 1:(d-1)) {
    
    reduced_nu = simulate_nu(eigs[1:j], d=j, alpha=alpha, a_sims=a_sims, nu_sims=nu_sims)
    
    # Cost-benefit calculation weighted by the eigenvalues of \Sigma
    benefit = sum(eigs[1:j] * (full_nu[1:j] - reduced_nu))
    cost = sum(eigs[(j+1):d] * (1 - full_nu[(j+1):d]))
    
    if ((benefit - cost) >= lb) {
      lb = benefit - cost
      weighted_k = j
    }
  }
  
  return(weighted_k)
}
