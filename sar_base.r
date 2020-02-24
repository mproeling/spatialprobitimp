require(Matrix)
require(spdep)

################################################################################
#
# Bayesian estimation of spatial autoregressive probit model (SAR probit)
# using MCMC sampling
#
# Stefan Wilhelm <Stefan.Wilhelm@financial.com>
# using code from Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
#
################################################################################

if (FALSE) {
  library(tmvtnorm)
  library(mvtnorm)
  library(Matrix)    # sparseMatrix
  source("sar_base.r")
  source("matrix_operations.r")
  source("stats_distributions.r")
  source("utility_functions.r")
}

# estimated tr(W^i) for i=1,...,100
# see LeSage (2009), chapter 4, p.97/98
#
# "The expectation of the quadratic form u'Au equals tr(A).
#  Since u_i^2 follows a chi^2 distribution with one degree of freedom."
# Pre-calculate traces tr(W^i) i=1,...,100 for the x-impacts calculations
#
# @param W spatial weight matrix (n x n)
# @param o highest order of W^i = W^o
# @param iiter number of MCMC iterations (we run this 50 times)
# @return (n x o) matrix with tr(W^i) in each column, for i=1..,o
tracesWi <- function(W, o=100, iiter=50) {
  n <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x o
  u  <- matrix(rnorm(n * iiter), nrow = n, ncol = iiter)   # (n x iiter)
  xx <- u
  trW_i[,1] <- apply(u * as.matrix(xx), 1, sum)    # tr(W^0) = 1
  for(i in 2:o ){
    xx <- W %*% xx  # (n x iter)
    trW_i[,i] <- apply(u * as.matrix(xx), 1, sum)  # u'(W^i)u; sum across all iterations
  }
  trW_i <- trW_i / iiter
  return(trW_i)
}

# Create mode function.
getmode <- function(x) {
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

# PURPOSE: draw rho from conditional distribution p(rho | beta, z, y)
# ---------------------------------------------------
#  USAGE:
#
#  results.rmin
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
#where epe0 <- t(e0) %*% e0
#      eped <- t(ed) %*% ed
#      epe0d<- t(ed) %*% e0
#
# detval1 umbenannt in => rho_grid = grid/vector of rhos
# detval2 umbenannt in => lndet = vector of log-determinants
# detval1sq = rvec^2
# lnbprior = vector of log prior densities for rho (default: Beta(1,1))
#
draw_rho <- function (rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d,
                      rho, nmk, nrho, lnbprior, u)
{
  # This is the scalar concentrated log-likelihood function value
  # lnL(rho). See eqn(3.7), p.48
  z <- epe0 - 2 * rho_grid * epe0d + rho_gridsq * eped
  z <- -nmk * log(z)
  den <- lndet + z + lnbprior          # vector of log posterior densities for rho vector
  # posterior density post(rho | data) \propto likelihood(data|rho,beta,z) * prior(rho)
  # log posterior density post(rho | data) \propto loglik(data|rho,beta,z) + log prior(rho)
  n <- nrho
  adj <- max(den)
  den <- den - adj                     # adjustieren/normieren der log density auf maximum 0;
  x <- exp(den)                        # density von p(rho) --> pdf
  isum <- sum(yy * (x[2:n] - x[1:(n - 1)])/2)
  z <- abs(x/isum)
  den <- cumsum(z)
  rnd <- u * sum(z)
  ind <- which(den <= rnd)
  idraw <- max(ind)
  if (idraw > 0 && idraw < nrho) {
    results <- rho_grid[idraw]
  }
  else {
    results <- rho
  }
  return(results)
}

beta_prior <- function(rvec,a1=1.01,a2=1.01){
  B          <- beta(a1,a2)
  num        <- (1+rvec)^(a1-1)
  num        <- num * (1-rvec)^(a2-1)
  den        <- 2^(a1+a2-1)
  results    <- (1/B)*num/den
  results[1] <- 0.001
  results[length(results)] <- 0.001
  return( results )
}

# Bayesian estimation of SAR probit model
#
# @param formula
sarprobit <- function(formula, W, data, subset, ...) {
  cl <- match.call()                     # cl ist object of class "call"
  mf <- match.call(expand.dots = FALSE)  # mf ist object of class "call"
  m  <- match(c("formula", "data", "subset"), names(mf), 0L)        # m is index vector
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())         # from here mf is a data.frame
  mt <- attr(mf, "terms")                # mt is object of class "terms" and "formula"
  y <- model.response(mf, "numeric")
  if (!is.null(W) && !is.numeric(W) && !inherits(W, "sparseMatrix") && nrow(W) != NROW(y))
    stop(gettextf("'W' must be a numeric square matrix, dimension %d should equal %d (number of observations)",
                  NROW(W), NROW(y)), domain = NA)

  X <- model.matrix(mt, mf, contrasts)
  sar_probit_mcmc(y, X, W, ...)
}

# faster update of matrix S = (I - rho * W) for new values of rho
#
# @param S template matrix of (I - rho * W)
# @param ind indizes to replaced
# @param W spatial weights matrix W
# @return (I - rho * W)
update_I_rW <- function(S, ind, rho, W) {
  S@x[ind] <- (-rho*W)@x
  return(S)
}

# PURPOSE: check if a matrix has an intercept in the
#          first row
# ---------------------------------------------------
#  USAGE: has_intercept( x )
#  where  x = n x K matrix
# ---------------------------------------------------
# RETURN:
#        TRUE - if the matrix has na intercept int the
#               first column
#        FALSE- if the matrix has no intercept
#        NULL - if the matrix has an intercept in a
#               column other than the first
has_intercept <- function( x ){
  n         <- nrow(x)
  k         <- ncol(x)
  intercept <- rep(1, n)
  results   <- FALSE
  for( c in 1:k ){
    check <- which( intercept == x[ ,c] )
    if( length( check ) == n && c == 1){
      results <- TRUE
    }else if( length( check ) == n && c != 1){
      print('has_intercept: intercept term must be in first column of the x-matrix')
      results <- NULL
    }
  }
  return( results )
}


sar_bayesian_estimation <- function(results){
  out       <- NULL
  #bayesian estimation
  bout_mean <- as.matrix(c(apply(results$bdraw,2,mean),mean(results$pdraw))) #parameter mean column
  bout_sd   <- as.matrix(c(apply(results$bdraw,2,sd)  ,sd(results$pdraw))) #parameter sd colum
  bout_sig  <- matrix(data=NA, nrow=nrow(bout_mean),ncol=1)
  #build bayesian significance levels
  draws     <- cbind( results$bdraw, results$pdraw )
  for( i in 1:ncol(draws) ){
    if( bout_mean[i,1] > 0){
      cnt <- which( draws[,i] > 0 )
    }else{
      cnt <- which( draws[,i] < 0 )
    }
    bout_sig[i,1] <- 1 - (length(cnt)/(results$ndraw-results$nomit))
  }

  out$bhat      <- bout_mean
  out$bhat_sd   <- bout_sd
  out$bhat_Bpval<- bout_sig
  out$bhat_names<- as.matrix( results$names )
  return(out)
}



#----------------------------------------------------
# PURPOSE: compute the eigenvalues for the weight matrix
# ---------------------------------------------------
#  USAGE: results = far_eigs(eflag, W)
#  results.rmin
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
# Adapted to R by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy
sar_eigs <- function(eflag, W){
  results <- NULL
  results$rmin <- -1
  results$rmax <-  1
  results$time <-  0
  # compute eigenvalues
  if( eflag == 1 ){
    t0           <- Sys.time()
    lambda       <- as.double( eigen( W, only.values=TRUE )$values )
    results$rmin <- 1/min(lambda)
    results$rmax <- 1
    results$time <- Sys.time() - t0
  }
  return( results )
}

#---------------------------------------------------
# PURPOSE: compute the log determinant |I_n - rho*W|
# using the user-selected (or default) method
# ---------------------------------------------------
# USAGE: detval = sar_lndet(lflag,W,rmin,rmax)
# where ldetflag,rmin,rmax,W contains input flags
# ---------------------------------------------------
# Code is now based on method spdep::do_ldet written by Roger Bivand
# ldetflag = 0 : Pace and Barry (1997) grid
# Pace, R. and Barry, R. (1997). Quick computation of spatial autoregressive estimators.
# Geographics Analysis, 29(3):232-247.
#
# ldetflag = 1 : Pace and LeSage (2004) Chebyshev approximation
# Pace, R. and LeSage, J. (2004). Chebyshev approximation of log-determinants of
# spatial weight matrices. Computational Statistics & Data Analysis, 45(2):179-
# 196.
#
# ldetflag = 2 : Barry and Pace (1999) MC approximation
# Barry, R. and Pace, R. (1999). Monte Carlo estimates of the log determinant of
# large sparse matrices. Linear Algebra and its Applications, 289(1-3):41-54.
sar_lndet <- function(ldetflag,W,rmin,rmax){
  results       <- NULL
  env <- new.env(parent=globalenv())
  listw <- mat2listw(W)
  assign("n", nrow(W), envir=env)
  assign("listw", listw, envir=env)
  assign("family", "SAR", envir=env)
  assign("similar", FALSE, envir=env)

  # do lndet approximation calculations if needed
  if( ldetflag == 0 ){ # no approximation
    t0            <- Sys.time()
    SE_classic_setup(env)
    # fine grid is already computed
    results$detval  <- get("detval1", env)
  }else if( ldetflag == 1 ) {
    t0   <- Sys.time()
    cheb_setup(env, q=4)  # Chebychev approximation, q = 4
    # For Chebyshev-Approximation we must compute the fine grid with do_ldet
    # to be later used for numerical integration
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else if( ldetflag == 2 ) {
    t0   <- Sys.time()
    mcdet_setup(env, p=16, m=30, which=1)
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else{
    # TODO: Pace and Barry, 1998 spline interpolation
    stop('sar_lndet: method not implemented')
  }
  results$time   <- Sys.time() - t0
  return( results )
}

# PURPOSE: computes Pace and Barry's grid for log det(I-rho*W) using sparse matrices
# -----------------------------------------------------------------------
# USAGE: out = lndetfull(W,lmin,lmax)
# where:
#             lmin  = lower bound on rho
#             lmax  = upper bound on rho
# -----------------------------------------------------------------------
# RETURNS: out = a structure variable
#          out.lndet = a vector of log determinants for 0 < rho < 1
#          out.rho   = a vector of rho values associated with lndet values
# -----------------------------------------------------------------------
# NOTES: should use 1/lambda(max) to 1/lambda(min) for all possible rho values
# -----------------------------------------------------------------------
# References: % R. Kelley Pace and  Ronald Barry. 1997. ``Quick
# Computation of Spatial Autoregressive Estimators'', Geographical Analysis
# -----------------------------------------------------------------------
#Written by:
#James P. LeSage, Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu
#
# Adapted to R by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy

# @param W spatial weight matrix
# @param rmin
# @param rmax
lndetfull <- function( W, rmin, rmax ){
  rvec    <- seq(rmin, rmax, 0.01)
  In      <- speye(nrow(W))
  niter   <- length(rvec)

  results       <- NULL
  results$rho   <- NULL
  results$lndet <- NULL

  for(i in 1:niter){
    rho             <- rvec[i]
    z               <- In - rho*W
    results$rho[i]  <- rho
    results$lndet[i]<- as.double(determinant(z, logarithm=TRUE)$modulus)
  }
  return(results)
}


# PURPOSE: This function computes the a
# the chebyshev a log determinant approximation.
# ---------------------------------------------------
#  USAGE:
#  where:
# ---------------------------------------------------
# SOURCE:
# Pace and LeSage 2004, Chebyshev Approximation of
# log-determinants of spatial weight matrices,
# Computational Statistics and Data Analysis, vol 45
#
# Written by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy

lndetChebyshev <- function(W, rmin=-1, rmax=1){
  results       <- NULL
  results$rho   <- NULL
  results$lndet <- NULL

  n  <- nrow(W)
  D1 <- W           #keep the notation of the paper
  D2 <- D1 %*% D1

  #pre-compute all required traces
  tr <- c(
    n,               #tr(I)
    0,               #tr(D)
    sum(D1 * t(D1)), #tr(D^2)
    sum(D1 * t(D2)), #tr(D^3)
    sum(D2 * t(D2))  #tr(D^4)
  )

  #pre-compute chebyshev approximation
  trTD <- c(
    tr[1],              #tr( TD0 ) = tr(I)
    tr[2],              #tr( TD1 ) = 0
    2*tr[3] - tr[1]  ,  #tr( TD2 ) = 2 * tr(TD2) - tr(I)
    4*tr[4] - 3*tr[2],  #tr( TD3 ) = 4 * tr(TD3) - 3 * tr(TD)
    8*tr[5] - 8*tr[3] + tr[ 1 ]
  )

  q    <- 4
  rvec <- seq(rmin, rmax, 0.01)
  niter<- length(rvec)
  for( i in 1:niter ){
    alpha <- rvec[ i ]
    c0    <- 0
    for( j in 1:(q+1) ){
      c0 <- c0 + chebyshev(j, q, alpha)*trTD[j]
    }
    results$rho[i]   <- alpha
    results$lndet[i] <- c0 - n/2 * chebyshev(1, q, alpha)
  }
  return( results )
}


# PURPOSE: This function is part of the computation of
# the chebyshev a log determinant approximation.
# ---------------------------------------------------
#  USAGE:
#  where:
# ---------------------------------------------------
# SOURCE:
# Pace and LeSage 2004, Chebyshev Approximation of
# log-determinants of spatial weight matrices,
# Computational Statistics and Data Analysis, vol 45
#----------------------------------------------------
# AUTHOR:
# Miguel Godinho de Matos
# Miguel.GodinhoMatos@gmail.com
# PhD Student of Engineering & Public Policy
# Carnegie Mellon University
chebyshev <- function(j, q, alpha){
  c0 <- 0
  for( k in 1:(q+1)){
    c0 <- c0 + log(1-(alpha*cos((pi*(k-0.5))/(q+1))))*cos((pi*(j-1)*(k-0.5))/(q+1))
  }
  c1 <- (2/(q+1)) * c0
  return(c1)
}

# % PURPOSE: computes and prints posterior model probabilities using log-marginals
# % ---------------------------------------------------
# %  USAGE: probs = model_probs(results1,results2,results3, ...)
# %  where: results_matrix is a  matrix of results structures returned by estimation
# %         functions, sar_g, sdm_g, sac_g, sem_g
# % e.g. result1 = sar_g(y,x,W,ndraw,nomit,prior);
# %      result2 = sem_g(y,x,W,ndraw,nomit,prior);
# % results_matrix = [result1 result2];
# % model_probs(results_matrix);
# % ---------------------------------------------------
# %  RETURNS: probs = a vector of posterior model probabilities
# % ---------------------------------------------------
#
# % written by:
# % James P. LeSage, 7/2003
# % Dept of Economics
# % University of Toledo
# % 2801 W. Bancroft St,
# % Toledo, OH 43606
# % jlesage@spatial-econometrics.com
# Ported to R by Miguel Godinho de Matos
#
# @param detval
# @param result_list a list of estimated model objects
# @return a vector of posterior model probabilities
model_probs <- function(detval, result_list){
  stopifnot(is.list(result_list))
  nmodels   <- length(result_list)
  nrho      <- nrow(detval)
  lmarginal <- matrix(ncol=nmodels, nrow=nrho, data=0)

  for(i in 1:nmodels){
    lmarginal[,i] <- result_list[[i]]$mlike
  }

  xx  <- exp(lmarginal- max(lmarginal))
  #yy  <- repmat(as.matrix(detval[,1]), nmodels)
  yy  <- matrix(detval[,1], length(detval[,1]), nmodels)

  isum<- colSums((yy[2:nrho,] + yy[1:(nrho-1),])*(xx[2:nrho,]-xx[1:(nrho-1),])/2)
  psum<- sum( isum )
  return( as.matrix( isum/psum ) )
}

sar_parse = function(prior, k, rho = 0.39, rhosd = 0.16){
  # PURPOSE: parses input arguments for sar_g models
  # ---------------------------------------------------
  #  USAGE: [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,c,T,inform_flag,a1,a2] =
  #                           sar_parse(prior,k)
  # where info contains the structure variable with inputs
  # and the outputs are either user-inputs or default values
  # ---------------------------------------------------

  # written by:
  # James P. LeSage, last updated 3/2010
  # Dept of Finance & Economics
  # Texas State University-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com


  # set defaults
  output = c()
  output$eflag = 0;     # default to not computing eigenvalues
  output$ldetflag = 1;  # default to 1999 Pace and Barry MC determinant approx
  output$mflag = 1;     # default to compute log marginal likelihood
  output$order = 50;    # there are parameters used by the MC det approx
  output$iter = 30;     # defaults based on Pace and Barry recommendation
  output$rmin = -1;     # use -1,1 rho interval as default
  output$rmax = 1;
  output$detval = 0;    # just a flag
  output$rho = rho;
  output$rhosd = rhosd;
  output$sige = 1.0;
  output$rval = 4;
  output$mm = 0;
  output$kk = 0;
  output$nu = 0;
  output$d0 = 0;
  output$a1 = 1.0;
  output$a2 = 1.0;
  output$diffuseprior = rep(0, k);   # diffuse prior for beta
  output$diagT = diag(k)*1e+12;
  output$prior_beta = 0;   # flag for diffuse prior on beta
  output$novi_flag = 0; # do vi-estimates
  output$inform_flag = 0;

  fields = names(prior)
  nf = length(fields)
  if(nf > 0){
    for(i in 1:nf){
      if (fields[i] == 'nu'){
        output$nu = prior$nu
      } else if(fields[i] == 'd0'){
        output$d0 = prior$d0
      } else if(fields[i] == 'rval'){
        output$rval = prior$rval
      } else if(fields[i] == 'a1'){
        output$a1 = prior$a1
      } else if(fields[i] == 'a2'){
        output$a2 = prior$a2
      } else if(fields[i] == 'm'){
        mm = prior$m
        kk = prior$k
        output$rval = gamm_rnd(1,1,mm,kk) # initial value for rval
      } else if(fields[i] == 'beta'){
        output$diffuseprior = prior$beta
        output$inform_flag = 1 # flag for informative prior on beta
      } else if(fields[i] == 'bcov'){
        output$diagT = prior$bcov
        output$inform_flag = 1
      } else if(fields[i] == 'rmin'){
        output$rmin = prior$rmin
        output$eflag = 0
      } else if(fields[i] == 'rmax'){
        output$rmax = prior$rmax
        output$eflag = 0
      } else if(fields[i] == 'lndet'){
        detval = prior$lndet
        output$ldetflag = -1
        output$eflag = 0
        output$rmin = detval[1]
        nr = length(detval)
        output$rmax = detval[nr]
      } else if(fields[i] == 'lflag'){
        tst = prior$lflag;
        if(tst == 0){
          output$ldetflag = 0
        } else if(tst == 1){
          output$ldetflag = 1
        } else if(tst == 2){
          output$ldetflag = 2
        } else {
          stop('sar_g: unrecognizable lflag value on input');
        }
      } else if(fields[i] == 'order'){
        output$order = prior$order
      } else if(fields[i] == 'iter'){
        output$iter = prior$iter
      } else if(fields[i] == 'novi'){
        output$novi_flag = prior$novi
      } else if(fields[i] == 'eig'){
        output$eflag = prior$eig
      }
    }
    return(output)
  } else {
    # the user has input a blank info structure
    # so we use the defaults
    return(output)
  }
}


sar_parse_multiple = function(prior, k, xobs, z, rho = 0.39, rhosd = 0.16){
  # PURPOSE: parses input arguments for sar_g models
  # ---------------------------------------------------
  #  USAGE: [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,c,T,inform_flag,a1,a2] =
  #                           sar_parse(prior,k)
  # where info contains the structure variable with inputs
  # and the outputs are either user-inputs or default values
  # ---------------------------------------------------

  # written by:
  # James P. LeSage, last updated 3/2010
  # Dept of Finance & Economics
  # Texas State University-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com
  # set defaults
  Y = z-1
  output = c()
  output$eflag = 0;     # default to not computing eigenvalues
  output$ldetflag = 1;  # default to 1999 Pace and Barry MC determinant approx
  output$mflag = 1;     # default to compute log marginal likelihood
  output$order = 50;    # there are parameters used by the MC det approx
  output$iter = 30;     # defaults based on Pace and Barry recommendation
  output$rmin = -1;     # use -1,1 rho interval as default
  output$rmax = 1;
  output$detval = 0;    # just a flag
  output$rho = rho;
  output$rhosd = rhosd;
  output$sige = 1.0;
  output$rval = 4;
  output$mm = 0;
  output$kk = 0;
  output$nu = 0;
  output$d0 = 0;
  output$a1 = 1.0;
  output$a2 = 1.0;
  output$diffuseprior = rep(0, Y+k);   # diffuse prior for beta
  output$diagT = diag(k+Y)*1e+12;
  output$prior_beta = 0;   # flag for diffuse prior on beta
  output$novi_flag = 0; # do vi-estimates
  output$inform_flag = 0;
  output$z = z;

  fields = names(prior)
  nf = length(fields)
  if(nf > 0){
    for(i in 1:nf){
      if (fields[i] == 'nu'){
        output$nu = prior$nu
      } else if(fields[i] == 'd0'){
        output$d0 = prior$d0
      } else if(fields[i] == 'rval'){
        output$rval = prior$rval
      } else if(fields[i] == 'a1'){
        output$a1 = prior$a1
      } else if(fields[i] == 'a2'){
        output$a2 = prior$a2
      } else if(fields[i] == 'm'){
        mm = prior$m
        kk = prior$k
        output$rval = gamm_rnd(1,1,mm,kk) # initial value for rval
      } else if(fields[i] == 'beta'){
        output$diffuseprior = c(rep(0,Y), prior$beta)
        output$inform_flag = 1 # flag for informative prior on beta
      } else if(fields[i] == 'bcov'){
        xobs = as.matrix(xobs)
        output$diagT = t(xobs)%*%xobs
        output$inform_flag = 1
      } else if(fields[i] == 'rmin'){
        output$rmin = prior$rmin
        output$eflag = 0
      } else if(fields[i] == 'rmax'){
        output$rmax = prior$rmax
        output$eflag = 0
      } else if(fields[i] == 'lndet'){
        detval = prior$lndet
        output$ldetflag = -1
        output$eflag = 0
        output$rmin = detval[1]
        nr = length(detval)
        output$rmax = detval[nr]
      } else if(fields[i] == 'lflag'){
        tst = prior$lflag;
          if(tst == 0){
          output$ldetflag = 0
        } else if(tst == 1){
          output$ldetflag = 1
        } else if(tst == 2){
          output$ldetflag = 2
        } else {
          stop('sar_g: unrecognizable lflag value on input');
        }
      } else if(fields[i] == 'order'){
        output$order = prior$order
      } else if(fields[i] == 'iter'){
        output$iter = prior$iter
      } else if(fields[i] == 'novi'){
        output$novi_flag = prior$novi
      } else if(fields[i] == 'eig'){
        output$eflag = prior$eig
      }
    }
    return(output)
  } else {
    # the user has input a blank info structure
    # so we use the defaults
    return(output)
  }
}

sar_marginal = function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2){
# PURPOSE: returns a vector of the log-marginal over a grid of rho-values
# -------------------------------------------------------------------------
  # USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
# where:       detval = an ngrid x 2 matrix with rho-values and lndet values
#                  e0 = y - x*b0;
#                 ed = Wy - x*bd;
#               epe0 = e0'*e0;
#               eped = ed'*ed;
#              epe0d = ed'*e0;
#               nobs = # of observations
#               nvar = # of explanatory variables
#            logdetx = log(det(x'*x))
#                 a1 = parameter for beta prior on rho
#                 a2 = parameter for beta prior on rho
# -------------------------------------------------------------------------
  # RETURNS: out = a structure variable
  #        out = log marginal, a vector the length of detval
  # -------------------------------------------------------------------------
    # NOTES: -this does not take any prior on beta, sigma into account, uses diffuse priors
  #         see sar_marginal2()
  # -------------------------------------------------------------------------

    # written by:
    # James P. LeSage, last updated 3/2010
  # Dept of Finance & Economics
  # Texas State University-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com
  if(is.matrix(detval)){n = nrow(detval)} else {n = length(detval)}
  nmk = (nobs-nvar)/2
  # C is a constant of integration that can vary with nvars, so for model
  # comparisions involving different nvars we need to include this
  bprior = beta_prior(detval[,1],a1,a2)
  bigC = log(bprior) + log(gamma(nmk)) - nmk*log(2*pi) - 0.5*logdetx
  iota = rep(1, n)
  z = as.numeric(epe0)*iota - 2*detval[,1]*as.numeric(epe0d) + detval[,1]*detval[,1]*as.numeric(eped)
  den = bigC + detval[,2] - nmk*log(z)
  out = Re(as.vector(den))
  return(out)
}

sar_marginal_multiple = function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2){
  # PURPOSE: returns a vector of the log-marginal over a grid of rho-values
  # -------------------------------------------------------------------------
  # USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
  # where:       detval = an ngrid x 2 matrix with rho-values and lndet values
  #                  e0 = y - x*b0;
  #                 ed = Wy - x*bd;
  #               epe0 = e0'*e0;
  #               eped = ed'*ed;
  #              epe0d = ed'*e0;
  #               nobs = # of observations
  #               nvar = # of explanatory variables
  #            logdetx = log(det(x'*x))
  #                 a1 = parameter for beta prior on rho
  #                 a2 = parameter for beta prior on rho
  # -------------------------------------------------------------------------
  # RETURNS: out = a structure variable
  #        out = log marginal, a vector the length of detval
  # -------------------------------------------------------------------------
  # NOTES: -this does not take any prior on beta, sigma into account, uses diffuse priors
  #         see sar_marginal2()
  # -------------------------------------------------------------------------

  # written by:
  # James P. LeSage, last updated 3/2010
  # Dept of Finance & Economics
  # Texas State University-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com
  if(is.matrix(detval)){n = nrow(detval)} else {n = length(detval)}
  nmk = (nobs-nvar)/2
  # C is a constant of integration that can vary with nvars, so for model
  # comparisions involving different nvars we need to include this
  bprior = beta_prior(detval[,1],a1,a2)
  bigC = log(bprior) + lgamma(nmk) - nmk*log(2*pi) - 0.5*logdetx
  iota = matrix(nrow=n, ncol=1, data=1)
  z = as.numeric(epe0)*iota - 2*detval[,1] * as.numeric(epe0d) + (detval[,1] * detval[,1]) * as.numeric(eped)
  den = bigC + detval[,2] - nmk*log(z)
  return(Re(as.vector(den)))
}

sar_marginal2 = function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,diffuseprior,TI,xs,ys,sige,W){
# PURPOSE: returns a vector of the log-marginal over a grid of rho-values
#          for the case of an informative prior on beta
# -------------------------------------------------------------------------
  # USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
# where:       detval = an ngrid x 2 matrix with rho-values and lndet values
#                  e0 = y - x*b0;
#                 ed = Wy - x*bd;
#               epe0 = e0'*e0;
#               eped = ed'*ed;
#              epe0d = ed'*e0;
#               nobs = # of observations
#               nvar = # of explanatory variables
#            logdetx = log(det(x'*x))
#                 a1 = parameter for beta prior on rho
#                 a2 = parameter for beta prior on rho
#                 c = prior mean for beta
#                TI = prior var-cov for beta
#                xs = x*sqrt(V) or x if homoscedastic model
#                ys = y*sqrt(V) or y is homoscedastic model
# -------------------------------------------------------------------------
  # RETURNS: out = a structure variable
  #        out = log marginal, a vector the length of detval
  # -------------------------------------------------------------------------
    # NOTES: -this is only an approximation based on the posterior mean Vi-estimates
  # -------------------------------------------------------------------------

    # written by:
    # James P. LeSage, last updated 3/2010
  # Dept of Finance & Economics
  # Texas State University-San Marcos
  # 601 University Drive
  # San Marcos, TX 78666
  # jlesage@spatial-econometrics.com

  if(is.matrix(detval)){n = nrow(detval)} else {n = length(detval)}
  nmk = (nobs-nvar)/2
  # C is a constant of integration that can vary with nvars, so for model
  # comparisions involving different nvars we need to include this
  bprior = beta_prior(detval[,1],a1,a2)
  bigC = log(bprior) + log(gamma(nmk)) - nmk*log(2*pi)
  iota = rep(1, n)
  z = as.numeric(epe0)*iota - 2*detval[,1]*as.numeric(epe0d) + detval[,1]*detval[,1]*as.numeric(eped)

  # add quadratic terms based on prior for beta
  Q1 = rep(0, n)
  Q2 = rep(0, n)
  xpxi = solve(t(xs)%*%xs)
  sTI = c(sige)*TI
  xpxis = solve(t(xs)%*%xs + sTI)
  logdetx = log(det(xpxis))
  bigC = bigC - 0.5*logdetx
  for(i in 1:n){
    rho = detval[i,1]
    D = Diagonal(nobs) - rho*W
    bhat = xpxi%*%(t(xs)%*%D%*%ys)
    beta = xpxis%*%(t(xs)%*%D%*%ys + sTI%*%diffuseprior)
    Q1[i] = t(diffuseprior - beta)%*%sTI%*%(diffuseprior - beta)
    Q2[i] = t(bhat - beta)%*%(t(xs)%*%xs)%*%(bhat - beta)
  }

  den = bigC + detval[,2] - nmk*log(z + Q1 + Q2)
  out = Re(as.vector(den))
  return(out)
}



marginal.effects <- function (object, ...)
  UseMethod("marginal.effects")

# compute marginal effects for every MCMC iteration of the
# estimated SAR probit model
marginal.effects.sarprobit <- function(object, o=100, ...) {
  # check for class "sarprobit"
  if (!inherits(object, "sarprobit"))
    stop("use only with \"sarprobit\" objects")

  nobs      <- object$nobs
  nvar      <- object$nvar   # number of explanatory variables
  p         <- ifelse(object$cflag == 0, nvar, nvar - 1) # number of non-constant variables
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  betadraws <- object$bdraw
  rhodraws  <- object$pdraw
  X         <- object$X # data matrix
  W         <- object$W # spatial weight matrix
  I_n       <- sparseMatrix(i=1:nobs, j=1:nobs, x=1) # sparse identity matrix

  # Monte Carlo estimation of tr(W^i) for i = 1..o, tr(W^i) is (n x o) matrix
  trW.i <- tracesWi(W, o=o, iiter=50)

  # Matrices (n x k) for average direct, indirect and total effects for all k (non-constant) explanatory variables
  # TODO: check for constant variables
  D <- matrix(NA, ndraw, p)
  I <- matrix(NA, ndraw, p)
  T <- matrix(NA, ndraw, p)

  # names of non-constant parameters
  if(object$cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }
  colnames(T) <- namesNonConstantParams
  colnames(D) <- namesNonConstantParams
  colnames(I) <- namesNonConstantParams

  ones <- rep(1, nobs)

  # loop all MCMC draws
  for(i in 1:ndraw) {

    # get parameters for this MCMC iteration
    # TODO: check for constant variables
    beta <- betadraws[i,]
    beff <- betadraws[i,-1]
    rho  <- rhodraws[i]

    # solving equation
    # (I_n - rho * W) mu = X beta
    # instead of inverting S = I_n - rho * W as in mu = (In -  rho W)^{-1} X beta.
    # QR-decomposition for sparse matrices
    S <- I_n - rho * W
    QR <- qr(S)  # class "sparseQR"
    mu <- qr.coef(QR, X %*% beta)      # n x 1

    # Marginal Effects for SAR Probit Model:
    # LeSage (2009), equation (10.10), p.294:
    # d E[y | x_r] / dx_r' = phi(S^{-1} I_n mean(x_r) beta_r) * S^{-1} I_n beta_r
    # This gives a (n x n) matrix
    # average direct effects = n^{-1} tr(S_r(W))

    # compute effects estimates (direct and indirect impacts) in each MCMC iteration
    rhovec <- rho^(0:(o-1)) # vector (o x 1) with [1, rho^1, rho^2 ..., rho^(o-1)],
    # see LeSage(2009), eqn (4.145), p.115

    # See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
    #   direct: M_r(D) = n^{-1} tr(S_r(W))           # efficient approaches available, see LeSage (2009), chapter 4, pp.114/115
    #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # S_r(W) is dense n x n matrix!
    # indirect: M_r(I) = M_r(T) - M_r(D)             # Computation of dense n x n matrix M_r(T) for total effects required!
    # See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
    pdfz <- dnorm(as.numeric(mu))                     # standard normal pdf phi(mu) = phi( (In -  rho W)^{-1} X beta )  # (n x 1)
    dd   <- sparseMatrix(i=1:nobs, j=1:nobs, x=pdfz)  # dd is diagonal matrix with pdfz as diagonal (n x n)

    dir      <- as.double(t(pdfz) %*% trW.i %*% rhovec /nobs)  # (1 x n) * (n x o) * (o x 1) = (1 x 1)
    # direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * (In -  rho W)^{-1} * beta_r
    avg_direct     <- dir * beff      # (p x 1)

    # We compute the average total effects without inverting S
    # as in the LeSage Matlab Code,
    # but using the QR decomposition of S which we already have!
    # average total effects = n^(-1) * 1_n' %*% (D %*% S^(-1) * b[r]) %*% 1_n
    #                       = n^(-1) * 1_n' %*% (D %*% x) * b[r]
    # where D=dd is the diagonal matrix containing phi(mu)
    # and x is the solution of S %*% x = 1_n, obtained from the QR-decompositon
    # of S. The average total effects is then the mean of (D %*% x) * b[r]
    # average total effects, which can be furthermore done for all b[r] in one operation.
    avg_total    <- mean(dd %*% qr.coef(QR, ones)) * beff
    avg_indirect       <- avg_total - avg_direct    # (p x 1)

    D[i,] <- avg_direct
    I[i,] <- avg_indirect
    T[i,] <- avg_total
  }

  summaryMarginalEffects <- function(x) {
    r <- cbind(
      apply(x, 2, mean),
      apply(x, 2, sd),
      apply(x, 2, mean)/apply(x, 2, sd))
    colnames(r) <- c("marginal.effect", "standard.error","z.ratio")
    return(r)
  }
  summary_direct <- summaryMarginalEffects(D)
  summary_indirect <- summaryMarginalEffects(I)
  summary_total <- summaryMarginalEffects(T)

  return(list(direct=D, indirect=I, total=T,
              summary_direct=summary_direct,
              summary_indirect=summary_indirect,
              summary_total=summary_total)
  )
}

# summary method for class "sarprobit"
summary.sarprobit <- function(object, var_names=NULL, file=NULL,
                              digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "sarprobit"))
    stop("use only with \"sarprobit\" objects")

  nobs      <- object$nobs
  nvar      <- object$nvar
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  draws     <- object$B

  #bayesian estimation
  bout_mean <- object$coefficients                         #parameter mean column
  bout_sd   <- apply(draws, 2, sd)                         #parameter sd colum
  # build bayesian significance levels
  # for parameter > 0 count all draws > 0  and determine P(X <= 0)
  # for parameter <= 0 count all draws <0  and determine P(X >= 0)
  bout_sig <- 1 - apply(draws, 2, function(x) { ifelse (mean(x) > 0, sum(x > 0), sum(x < 0)) }) / ndraw
  #standard asymptotic measures
  bout_t    <- bout_mean / bout_sd             #t score b/se
  bout_tPval<- (1 - pt( abs(bout_t), nobs ))*2 #two tailed test = zero probability = z-prob
  #name definition
  if( is.null(var_names)){
    bout_names<- as.matrix(object$names)
  }else{
    bout_names<- as.matrix(var_names)
  }

  if(is.null(file)){file <- ""}#output to the console
  #HEADER
  write(sprintf("--------MCMC spatial autoregressive probit--------"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f %s", object$time, attr(object$time, "units"))  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("# of 0 Y values = %6d, # of 1 Y values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names,
                                 c("Estimate", "Std. Dev", "p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
               signif.stars = getOption("show.signif.stars"))
  #if (getOption("show.signif.stars")) {
  # The problem: The significance code legend printed by printCoefMat()
  # is too wide for two-column layout and it will not wrap lines...
  # The solution: using cat() instead of print() and use line breaks
  # cat(paste(strwrap(x, width = 70), collapse = "\\\\\n"), "\n")
  # http://r.789695.n4.nabble.com/Sweave-line-breaks-td2307755.html
  #  Signif <- symnum(1e-6, corr = FALSE, na = FALSE,
  #              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  #              symbols = c("***", "**", "*", ".", " "))
  #x <- paste("Signif. codes: ", attr(Signif, "legend"), "\n", sep="")
  #cat(paste(strwrap(x, width = getOption("width")), collapse = "\\\n"), "\n")
  #}
  return(invisible(coefficients))
}

# prints the marginal effects/impacts of the SAR probit fit
impacts.sarprobit <- function(obj, file=NULL,
                              digits = max(3, getOption("digits")-3), ...) {
  if (!inherits(obj, "sarprobit"))
    stop("use only with \"sarprobit\" objects")

  if(is.null(file)){file <- ""}#output to the console
  write(sprintf("--------Marginal Effects--------"), file, append=T)
  write(sprintf(""), file, append=T)

  #(a) Direct effects
  write(sprintf("(a) Direct effects"), file, append=T)
  direct <- cbind(
    lower_005=apply(obj$direct, 2, quantile, prob=0.05),
    posterior_mean=colMeans(obj$direct),
    upper_095=apply(obj$direct, 2, quantile, prob=0.95)
  )
  printCoefmat(direct, digits = digits)

  # (b) Indirect effects
  write(sprintf(""), file, append=T)
  write(sprintf("(b) Indirect effects"), file, append=T)
  indirect <- cbind(
    lower_005=apply(obj$indirect, 2, quantile, prob=0.05),
    posterior_mean=colMeans(obj$indirect),
    upper_095=apply(obj$indirect, 2, quantile, prob=0.95)
  )
  printCoefmat(indirect, digits = digits)

  # (c) Total effects
  write(sprintf(""), file, append=T)
  write(sprintf("(c) Total effects"), file, append=T)
  total <- cbind(
    lower_005=apply(obj$total, 2, quantile, prob=0.05),
    posterior_mean=colMeans(obj$total),
    upper_095=apply(obj$total, 2, quantile, prob=0.95)
  )
  printCoefmat(total, digits = digits)

  return(invisible(list(direct=direct, indirect=indirect, total=total)))
}

# concatenate two objects of class "sarprobit"
# c.sarprobit works in the same way as boot:::c.boot().
c.sarprobit <- function(...) {
  args <- list(...)
  nm <- lapply(args, names)
  if (!all(sapply(nm, function(x) identical(x, nm[[1]]))))
    stop("arguments are not all the same type of \"sarprobit\" object")
  res <- args[[1]]
  res$time  <- as.difftime(max(as.numeric(sapply(args, "[[", "time"), units="secs")), units="secs")
  res$ndraw <- sum(sapply(args, "[[", "ndraw"))
  res$nomit <- sum(sapply(args, "[[", "nomit"))
  res$B     <- do.call(rbind, lapply(args, "[[", "B"))
  res$bdraw <- do.call(rbind, lapply(args, "[[", "bdraw"))
  res$pdraw <- do.call(rbind, lapply(args, "[[", "pdraw"))
  res$coefficients <- colMeans(res$B)
  res$beta  <- res$coefficients[1:res$nvar]
  res$rho   <- res$coefficients[res$nvar+1]
  res
}


# extract the coefficients
coef.sarprobit <- function(object, ...) {
  if (!inherits(object, "sarprobit"))
    stop("use only with \"sarprobit\" objects")
  return(object$coefficients)
}

# extract the coefficients
coefficients.sarprobit <- function(object, ...) {
  UseMethod("coef", object)
}

# plot MCMC results for class "sarprobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.sarprobit <- function(x, which=c(1, 2, 3),
                           ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
  if (!inherits(x, "sarprobit"))
    stop("use only with \"sarprobit\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 3))
    stop("'which' must be in 1:3")

  names <- x$names
  B <- x$B
  k <- ncol(B)

  show <- rep(FALSE, 3)
  show[which] <- TRUE

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (show[1L]) {
    # trace plots
    for (i in 1:k) {
      plot(1:nrow(B), B[,i], type="l", xlab="iteration", ylab=names[i], main=substitute("Trace plot of "*x, list(x=names[i])), ...)
      if (!is.null(trueparam)) abline(h=trueparam[i], col="red", lty=2)
    }
  }

  if (show[2L]) {
    # ACFs
    for (i in 1:k) {
      acf(B[,i], main=substitute("ACF of "*x, list(x=names[i])), ...)
    }
  }

  if (show[3L]) {
    # posterior distribution
    for (i in 1:k) {
      plot(density(B[,i]), main=substitute("Posterior distribution of "*x, list(x=names[i])), ...)
      if (!is.null(trueparam)) abline(v=trueparam[i], col="red", lty=2)
    }
  }
}

# return fitted values of SAR probit (on response scale vs. linear predictor scale)
# TODO: see fitted.lm() for comparison
fitted.sarprobit <- function(object, ...) {
  object$fitted.value
}

# Extract Log-Likelihood; see logLik.glm() for comparison
# Method returns object of class "logLik" with at least one attribute "df"
# giving the number of (estimated) parameters in the model.
# see Marsh (2000) equation (2.8), p.27
logLik.sarprobit <- function(object, ...) {
  X <- object$X
  y <- object$y
  n <- nrow(X)
  k <- ncol(X)
  W <- object$W
  beta <- object$beta
  rho <- object$rho
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
  S <- I_n - rho * W
  D <- diag(1/sqrt(diag(S %*% t(S))))  # D = diag(E[u u'])^{1/2}  (n x n)
  Xs <- D %*% solve(qr(S), X)          # X^{*} = D * (I_n - rho * W)^{-1} * X
  #Xs <- D %*% solve(S) %*% X
  F <- pnorm(as.double(Xs %*% beta))    # F(X^{*} beta)  # (n x 1)
  lnL <- sum(log(F[y == 1])) + sum(log((1 - F[y == 0]))) # see Marsh (2000), equation (2.8)
  #lnL <- sum(ifelse(y == 1, log(pnorm(xb)), log(1 - pnorm(xb))))
  out <- lnL
  class(out) <- "logLik"
  attr(out,"df") <- k+1                 # k parameters in beta, rho
  return(out)
}
