
sar_continuous_mcmc <- function(y, x, W, ndraw=1000, burn.in=100, thinning=1,
                                prior=list(a1=1, a2=1, c=rep(0, ncol(X)),
                                           T=diag(ncol(X))*1e12, lflag = 0),
                                start=list(rho=0.75, beta=rep(0, ncol(X))),
                                m=10, computeMarginalEffects=TRUE, showProgress=TRUE){

  x = X
  ndraw=2500
  burn.in=500
  thinning=1
  prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0)
  m=10
  computeMarginalEffects=TRUE
  showProgress=TRUE

  n  <- nrow(x)            # number of observations
  n1 <- nrow(x)
  k <- ncol(x)             # number of of parameters/exogenous variables
  n2 <- dim(W)[1]
  n4 <- dim(W)[2]

  if(!is.matrix(x)){x = as.matrix(x)}

  results = c()
  results$nobs  = n
  results$nvar  = k
  results$y = y

  if(nargs() == 5){
    prior$lflag = 1
  }

  priors = sar_parse(prior,k);

  # check if the user handled the intercept term okay
  n = length(y)
  if(sum(x[,1]) != n){
    tst = colSums(x)          # we may have no intercept term
    ind = tst[tst == n]       # we do have an intercept term
    if(length(ind) > 0){
      stop('sar_g: intercept term must be in first column of the x-matrix');
    } else if(length(ind) == 0){# case of no intercept term
      cflag = 0
      p = dim(x)[2]
    }
  } else if(sum(x[,1]) == n){ # we have an intercept in the right place
    cflag = 1
    p = dim(x)[2]-1
  }

  results$cflag = cflag
  results$p = p

  if(n1 != n2){
    stop('sar_g: wrong size weight matrix W');
  } else if(n1 != n){
    stop('sar_g: wrong size weight matrix W');
  }

  if(is.vector(y)){nchk = length(y)} else {nchk = dim(y)[1]}

  if(nchk != n){
    stop('sar_g: wrong size y vector input');
  }

  results$order = priors$order
  results$iter = priors$iter

  output1 =  sar_eigs(priors$eflag, W)
  rmin = output1$rmin
  rmax = output1$rmax
  results$time1 = output1$time

  output2 = sar_lndet(priors$ldetflag,W,rmin,rmax)
  detval = output2$detval
  results$time2 = output2$time

  results$order = priors$order
  results$iter = priors$iter

  # storage for draws
  bsave = matrix(0, ndraw-burn.in, k)
  if(priors$mm != 0){
    rsave = rep(0, ndraw-burn.in)
  }
  psave = rep(0, ndraw-burn.in)
  ssave = rep(0, ndraw-burn.in)
  vmean = rep(0, n)

  # ====== initializations
  # compute this stuff once to save time
  TI = solve(priors$diagT)
  TIc = TI%*%priors$diffuseprior

  in.ones = rep(1, n)
  V = in.ones
  vi = in.ones
  Wy = W%*%y  # W has been ordered so that the outcome Wy is also ordered (1 to 49)
  # in the original dataset wmat, the edgedata is not sorted so that results
  # in a sparse matrix with strange W orientation.

  # Some precalculated quantities for drawing rho from spatialprobit package
  # rho ~ Beta(a1, a2) prior
  rho_grid <- detval[,1]  # rho grid values
  lndet    <- detval[,2]  # log-determinant grid values
  lnbprior <- log(beta_prior(detval[,1], priors$a1, priors$a2))
  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval)  # do_ldet() gives only 2000 statt 2001 Gridpoints
  nmk      <- (n-k)/2
  rho_gridsq <- rho_grid * rho_grid
  yy       <- (rho_grid[2:nrho] + rho_grid[1:(nrho-1)])

  if(prior2$novi == 0){  # fit heteroscedastic model

    if(showProgress){pb <- txtProgressBar(min=0, max=(thinning * ndraw), initial=0, style=3)}
    iter = 1
    while(iter <= ndraw){ # start sampling;
      # update beta
      xs = as.matrix(x)*sqrt(V)
      ys = sqrt(V)*y
      Wys = sqrt(V)*Wy
      AI = qr.solve(t(xs)%*%xs + priors$sige*TI, diag(k))
      yss = ys - priors$rho*Wys
      xpy = t(xs)%*%yss
      b = t(xs)%*%yss + priors$sige*TIc
      b0 = qr.solve(t(xs)%*%xs + priors$sige*TI, b)
      bhat = norm_rnd(priors$sige*AI) + b0
      xb = xs%*%bhat

      # update sige
      nu1 = n + 2*priors$nu
      e = (yss - xb)
      d1 = 2*priors$d0 + t(e)%*%e
      chi =  rchisq(1,nu1)
      priors$sige = as.numeric(d1/chi)

      # update vi
      ev = y - priors$rho*Wy - as.matrix(x)%*%bhat
      chiv = rchisq(n,priors$rval+1)
      # chiv = chi2rnd(rval+1,n,1); % Statistics Toolbox function
      vi = ((ev*ev/priors$sige) + in.ones*priors$rval) / chiv
      V = in.ones/vi

      # update rval
      if(priors$mm != 0){
        priors$rval = gamm_rnd(1,1, priors$mm, priors$kk)
      }

      # we use griddy Gibbs to perform rho-draw
      b0 = qr.solve((t(xs)%*%xs + priors$sige*TI), (t(xs)%*%ys + priors$sige*TIc))
      bd = qr.solve((t(xs)%*%xs + priors$sige*TI), (t(xs)%*%Wys + priors$sige*TIc))
      e0 = ys - xs%*%b0
      ed = Wys - xs%*%bd
      epe0 = as.numeric(t(e0)%*%e0)
      eped = as.numeric(t(ed)%*%ed)
      epe0d = as.numeric(t(ed)%*%e0)
      priors$rho <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, priors$rho, nmk=nmk, nrho=nrho, lnbprior, u=u[iter + burn.in])

      if(iter > burn.in){ # if we are past burn-in, save the draws
        bsave[iter-burn.in,1:k] = t(bhat)
        ssave[iter-burn.in] = priors$sige
        psave[iter-burn.in] = priors$rho
        vmean = vmean + vi
        if(priors$mm != 0){
          rsave[iter-burn.in] = priors$rval
        }
      }

      iter = iter + 1
      if(showProgress){setTxtProgressBar(pb, iter)}
    } # End the sampling
    if(showProgress){close(pb)}
  } else if(prior2$novi == 1){  # fit homoscedastic model

    if(showProgress){pb <- txtProgressBar(min=0, max=(thinning * ndraw), initial=0, style=3)}
    iter = 1
    if(!is.matrix(x)){x=as.matrix(x)}
    xpx = t(x)%*%x
    xpy = t(x)%*%y
    Wy = W%*%y
    xpWy = t(x)%*%Wy


    while(iter <= ndraw){ # start sampling;
      AI = qr.solve(xpx + priors$sige*TI, diag(k))
      ys = y - priors$rho*Wy
      b = t(as.matrix(x))%*%ys + priors$sige*TIc
      b0 = qr.solve(xpx + priors$sige*TI, b)
      bhat = norm_rnd(priors$sige*AI) + b0
      xb = x%*%bhat

      # update sige
      nu1 = n + 2*priors$nu
      e = (ys - xb)
      d1 = 2*priors$d0 + t(e)%*%e
      chi =  rchisq(1,nu1)
      priors$sige = as.numeric(d1/chi)

      # update rho using griddy Gibbs
      AI = qr.solve(xpx + priors$sige*TI, diag(k))
      b0 = qr.solve(xpx + priors$sige*TI, xpy + priors$sige*TIc)
      bd = qr.solve(xpx + priors$sige*TI, xpWy + priors$sige*TIc)
      e0 = y - as.matrix(x)%*%b0
      ed = Wy - as.matrix(x)%*%bd
      epe0 = as.numeric(t(e0)%*%e0)
      eped = as.numeric(t(ed)%*%ed)
      epe0d = as.numeric(t(ed)%*%e0)
      priors$rho  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, priors$rho, nmk=nmk, nrho=nrho, lnbprior, u=u[iter + burn.in])

      if(iter > burn.in){ # if we are past burn-in, save the draws
        bsave[iter-burn.in,1:k] = t(bhat)
        ssave[iter-burn.in] = priors$sige
        psave[iter-burn.in] = priors$rho
        vmean = vmean + vi
      }

      iter = iter + 1
      if(showProgress){setTxtProgressBar(pb, iter)}
    } # End the sampling
    if(showProgress){close(pb)}
  } else {
    stop('unrecognized prior2.novi value on input,
         cannot decide whether to fit a homo- or heteroscedastic model')
  }

  # pre-calculate traces for the x-impacts calculations
  uiter = 50
  maxorderu = 100
  nobs = n
  rv = matrix(rnorm(nobs * uiter), nobs, uiter)
  tracew = rep(0, maxorderu)
  wjjju = rv
  for(jjj in 1:maxorderu){
    wjjju = W%*%wjjju
    tracew[jjj] = mean(mean(rv*wjjju))
  }

  traces = tracew
  traces[1] = 0
  traces[2] = sum(sum(t(W)*W))/nobs
  trs = c(1, traces)
  ntrs = length(trs)
  trbig = t(trs)

  if(cflag == 1){
    bdraws = bsave[,c(2:ncol(bsave))]
  } else if(cflag == 0){
    bdraws = bsave
  }
  pdraws = psave

  ree = seq(0, ntrs-1)

  rmat = rep(0, ntrs)
  total = array(rep(matrix(ndraw-burn.in*p, ndraw-burn.in, p), ntrs), dim = c(ndraw-burn.in,p,ntrs))
  direct = array(rep(matrix(ndraw-burn.in*p, ndraw-burn.in, p), ntrs), dim = c(ndraw-burn.in,p,ntrs))
  indirect = array(rep(matrix(ndraw-burn.in*p, ndraw-burn.in, p), ntrs), dim = c(ndraw-burn.in,p,ntrs))

  for(i in 1:(ndraw-burn.in)){
    rmat = pdraws[i]^ree
    for(j in 1:p){
      beta = bdraws[i,j]
      brmat = beta*rmat
      btrmat = (beta*trbig)*rmat
      for(k in 1:ntrs){
        total[i,j,k] = brmat[k]
        direct[i,j,k] = btrmat[k]
        indirect[i,j,k] = total[i,j,k] - direct[i,j,k]
      }
    }
  }

  # compute posterior means and log marginal likelihood for return arguments
  bmean = colMeans(bsave)
  beta = t(bmean)
  rho = mean(psave)
  sige = mean(ssave)
  vmean = vmean/(ndraw-burn.in)
  V = in.ones/vmean

  results$sige = sige
  nobs = dim(x)[1]
  nvar = dim(x)[2]
  xs = as.matrix(x)*sqrt(V)
  ys = sqrt(V)*y
  Wys = W%*%ys
  AI = solve(t(xs)%*%xs + sige*TI)
  b0 = AI%*%(t(xs)%*%ys + sige*TIc)
  bd = AI%*%(t(xs)%*%Wys + sige*TIc)
  e0 = ys - xs%*%b0
  ed = Wys - xs%*%bd
  epe0 = t(e0)%*%e0
  eped = t(ed)%*%ed
  epe0d = t(ed)%*%e0
  logdetx = log(det(t(xs)%*%xs + sige*TI))
  if(priors$inform_flag == 0){
    mlike = sar_marginal_multiple(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,priors$a1,priors$a2)
  } else if(priors$inform_flag == 1){
    mlike = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,priors$a1,priors$a2,priors$diffuseprior,TI,xs,ys,sige,W);
  }
  yhat = qr.solve(Diagonal(nobs) - rho*W, as.matrix(x)%*%t(beta))
  e = y - yhat

  # compute R-squared
  epe = t(e)%*%e
  sige = epe/(n-k)
  results$sigma = sige
  ym = y - mean(y)
  rsqr1 = epe
  rsqr2 = t(ym)%*%ym
  results$rsqr = 1- rsqr1/rsqr2     # r-squared
  rsqr1 = rsqr1/(nobs-nvar)
  rsqr2 = rsqr2/(nobs-1.0)
  results$rbar = 1 - (rsqr1/rsqr2)  # rbar-squared

  # write output
  results$meth  = 'sar_g'
  results$total = total
  results$direct = direct
  results$indirect = indirect
  results$beta_std = apply(bsave, 2, sd)
  results$sige_std = sd(ssave)
  results$rho_std = sd(psave)
  results$beta = beta
  results$rho = rho
  results$bdraw = bsave
  results$pdraw = psave
  results$sdraw = ssave
  results$mlike = mlike
  results$vmean = vmean
  results$yhat  = yhat
  results$resid = e
  results$bmean = priors$diffuseprior
  results$bstd  = sqrt(diag(priors$diagT))
  results$ndraw = ndraw
  results$burn.in = burn.in
  results$nu = priors$nu
  results$d0 = priors$d0
  results$a1 = priors$a1
  results$a2 = priors$a2
  results$tflag = 'plevel'
  results$rmax = priors$rmax
  results$rmin = priors$rmin
  results$lflag = priors$ldetflag
  results$lndet = detval
  results$novi  = prior2$novi
  results$priorb = priors$inform_flag
  results$W = W
  results$x = x
  results$y = y

  if(priors$mm != 0){
    results$rdraw = rsave
    results$m     = priors$mm
    results$k     = priors$kk
  }else{
    results$r     = priors$rval
    results$rdraw = 0
  }

  class(results)    <- "sarprobit"
  return(results)
}

