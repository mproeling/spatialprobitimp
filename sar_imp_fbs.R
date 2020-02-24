### this function is an adaption from sar_combined_mcmc from the spatialprobit package:
# https://github.com/cran/spatialprobit

##### Full Bayes MODEL #####

# see sar_imp_cut.R for an explanation

sar_combined_mcmc_imp_fbs = function(y, x, W, ndraw=1000, burn.in=100, thinning=1,
                               prior=list(a1=1, a2=1, beta = as.list(beta), T=diag(ncol(X))*1e12, lflag=0),
                               start=list(rho=0.39, rhosd=0.16, beta=rep(0, ncol(X)), phi=c(-Inf, 0:(max(y)-1), Inf)),
                               m=10, computeMarginalEffects=TRUE, showProgress=FALSE,
                               method=vector("character", meth.length = ncol(x)), model=1){

  # hier zit eigenlijk een bug in de definitie van beta in de functie, as list maakt hier een lijst van 4 delen met een 0, ipv 1 deel met 4*0.
  # dat gaat nu goed omdat beta al is gedefinieerd voor het runnen, maar is niet optimaal.
  #y = yf1
  #y = log(yf10)
  #x = as.matrix(X);  ndraw=1000;  burn.in=100;  thinning=1;  prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0);  rho = 0.39;
  #rhosd = 0.16;  m=10;  computeMarginalEffects=TRUE;  showProgress=TRUE;  method = c("continuous");
  #start=list(rho=0.39, rhosd=0.16, beta=rep(0, ncol(X))) #, phi=c(-Inf, 0:(max(y)-1), Inf))

  if(!is.data.frame(y)){y = as.data.frame(y)}
  if(!is.matrix(x)){x = as.matrix(x)}

  Yk = ncol(y)
  if(Yk != length(method)){
    stop("definition of method needs to occur for every covariate, so if Y1 is continuous and Y2 is dichotomous; method = c('continuous', 'probit')")
  }

  obs = list(); Nobs = list(); Nmiss = list(); vmean = list();
  in.ones = list(); nmk = list(); missing = list(); y_imp = list();
  S = list(); H = list(); mu = list(); trW.i = list()
  betadraws = list(); ind = list(); ind2 = list(); zmean = list();
  beta = list(); for(i in 1:ncol(y)){beta[[i]] = rep(0, (ncol(X) + (i-1)))}
  total = list(); direct = list(); indirect = list(); lower = list(); upper = list()

  for(i in 1:Yk){
    obs[[i]] = which(!is.na(y[,i]))
    missing[[i]] = which(is.na(y[,i]))
    Nobs[[i]] = length(which(!is.na(y[,i])))
    Nmiss[[i]] = length(which(is.na(y[,i])))
    y_imp[[i]] = matrix(0, ndraw*thinning, Nmiss[[i]])

    # replace the missing values in Y with random values for the first run so that
    # the parameter estimates are based on all observatoins and the imputed values
    # influence the paremeter estimates
    if(method[i] == "probit"){
      imp = runif(Nmiss[[i]])
      imp = ifelse(imp > 0.5, 1, 0)
      y[missing[[i]], i] <- imp
      rho = start$rho
    } else if (method[i] == "orderedprobit"){
      imp = rnorm(Nmiss[[i]], mean(y[[i]], na.rm = TRUE), sd(y[[i]], na.rm = TRUE))
      imp <- cut(as.double(imp), breaks=start$phi, labels=FALSE, ordered_result = TRUE) # split according to phi values
      y[missing[[i]], i] <- imp
      rho = start$rho
    } else {
      imp = rnorm(Nmiss[[i]], mean(y[[i]], na.rm = TRUE), sd(y[[i]], na.rm = TRUE))
      y[missing[[i]], i] <- imp
      rho = start$rho
    }

    # partition data
    assign(paste0("W", i), W) # select the full W instead of the observed only
    assign(paste0("y", i), y[,i]) # select the full y instead of the observed only
    assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i)))
    # W has been ordered so that the outcome Wy is also ordered (1 to 49)
    # in the original dataset wmat, the edgedata is not sorted so that results
    # in a sparse matrix with strange W orientation.

    n  <- nrow(x)  # number of observations
    n1 <- nrow(x)
    k <- ncol(x)              # number of of parameters/exogenous variables
    n2 <- dim(get(paste0("W", i)))[1]
    n4 <- dim(get(paste0("W", i)))[2]

    vmean[[i]] = rep(0, n) # this is based on n_obs so this can be different on every i
    in.ones[[i]] = rep(1, n) # this is based on n_obs so this can be different on every i
    nmk[[i]] = (n-k)/2

    nchk = list()
    if(is.vector(get(paste0("y", i)))){nchk[[i]] = length(get(paste0("y", i)))} else {nchk[[i]] = dim(get(paste0("y", i)))[1]}
    if(nchk[[i]] != n){
      stop('sar_g: wrong size y vector input');
    }

    x_obs = as.matrix(cbind(y[, c(0:(i-1))], x)) # add y

    if(i == 1){
      # create x_mis for the imputation, only necessary for Y1 with least miss because the other data for Y2-Yk are already in X_obs
      x_mis = x[missing[[i]],]
      # check if the user handled the intercept term okay
      # n = length(y)
      ind[[i]] <- match(nrow(X), apply(X,2,sum))
      if( is.na(ind[[i]]) ){
        stop('sarprobit: intercept term must be in first column of the X-matrix')
      } else if( ind[[i]] == 1 ){
        cflag <- 1
        p     <- k - 1
      }
    } else {
      colnames(x_obs)[1:i-1] = paste0("Y", c(1:(i-1)))
      x_obs = x_obs[,c(i, 1:(i-1), (i+1):ncol(x_obs))]
    }
    location.intercept = which(colnames(x_obs) == "Intercept")

    if(i > 1 & any(diff(as.numeric(Nobs)) > 0)){ # check whether the order of missingness in Y is correct
      stop('Y input seems incorrectly ordered. Order should be from least to most missingness from left to right (Y1_miss < Y2_miss < ... < Yk_miss)')
    }

    if(method[i] == "probit"){
      if(length(c(which(y[[i]] == 0 ), which(y[[i]] == 1))) != length(y[[i]])){
        stop('sarprobit: not all y-values are 0 or 1')
      }

      if (is.numeric(prior$c) && length(prior$c) == k) {
        c_prior <- prior$c
      }

      I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
      # prepare computation of (I_n - rho * W)
      if (class(get(paste0("W", i))) == "dgCMatrix") {
        I <- sparseMatrix(i=1:n, j=1:n, x=Inf)
        S[[i]] <- (I - rho * get(paste0("W", i)))
        ind[[i]]  <- which(is.infinite(S[[i]]@x))
        ind2[[i]] <- which(!is.infinite(S[[i]]@x))
        S[[i]]@x[ind[[i]]] <- 1
      } else {
        S[[i]] <- I_n - rho * get(paste0("W", i))
      }

      H[[i]] <- t(S[[i]]) %*% S[[i]]            # precision matrix H for beta | rho, z, y
      QR <- qr(S[[i]])                          # class "sparseQR"
      mu[[i]] <- solve(QR, x_obs %*% beta[[i]]) # this beta is the starting value

      # truncation points for z, depend only on y, can be precalculated
      lower[[i]] <- ifelse(get(paste0("y", i)) > 0, 0,  -Inf)
      upper[[i]] <- ifelse(get(paste0("y", i)) > 0, Inf,   0)

      assign(paste0("priors",i), sar_parse_multiple(prior, k, x_obs, i))

      xpx  <- t(x_obs) %*% x_obs       # (X'X)            # k x k
      xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
      xxpxI <- x_obs %*% xpxI          # X(X'X)^(-1)     # n x k (better, compromise)
      AA    <- solve(xpx + solve(get(paste0("priors", i))$diagT)) # (X'X + T^{-1})^{-1}

      # draw from multivariate normal beta ~ N(c, T). we can precalculate
      # betadraws ~ N(0, T) befor running the chain and later just create beta as
      # beta = c + betadraws ~ N(c, T)
      betadraws[[i]] <- rmvnorm(n=(burn.in + ndraw * thinning), mean=rep(0, k + (i-1)), sigma=AA)

      # matrices for direct and indirect impacts
      direct[[i]]   <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      indirect[[i]] <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      total[[i]]    <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      zmean[[i]]    <- matrix(0, n, 1)

      # names of non-constant parameters
      if(cflag == 0) {
        namesNonConstantParams <- colnames(x_obs)
      } else {
        namesNonConstantParams <- colnames(x_obs)[-location.intercept]
      }
      colnames(total[[i]]) <- namesNonConstantParams
      colnames(direct[[i]])   <- namesNonConstantParams
      colnames(indirect[[i]]) <- namesNonConstantParams

      if (computeMarginalEffects) {
        # simulate Monte Carlo estimation of tr(W^i) for i = 1..o before MCMC iterations
        trW.i[[i]] <- tracesWi(W, o=100, iiter=50)
      }
    } else if(method[i] == "orderedprobit"){
      J <- max(get(paste0("y", i)))
      if( sum(get(paste0("y", i)) %in% 1:J) != length( get(paste0("y", i)))  ){
        stop('sarorderedprobit: not all y-values are in 1...J')
      }

      c_prior <- rep(0, k)
      if (is.numeric(prior$c) && length(prior$c) == k) {
        c_prior <- prior$c
      }

      phi  <- start$phi
      if (is.null(phi)) {
        phi <- c(-Inf, 0:(J-1), +Inf)
      }

      I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
      # prepare computation of (I_n - rho * W)
      if (class(get(paste0("W", i))) == "dgCMatrix") {
        I <- sparseMatrix(i=1:n, j=1:n, x=Inf)
        S[[i]] <- (I - rho * get(paste0("W", i)))
        ind[[i]]  <- which(is.infinite(S[[i]]@x))
        ind2[[i]] <- which(!is.infinite(S[[i]]@x))
        S[[i]]@x[ind[[i]]] <- 1
      } else {
        S[[i]] <- I_n - rho * get(paste0("W", i))
      }

      H[[i]] <- t(S[[i]]) %*% S[[i]]            # precision matrix H for beta | rho, z, y
      QR <- qr(S[[i]])                          # class "sparseQR"
      mu[[i]] <- solve(QR, x_obs %*% beta[[i]]) # this beta is the starting value

      # truncation points for z, depend only on y, can be precalculated
      lower[[i]] <- ifelse(get(paste0("y", i)) > 0, 0,  -Inf)
      upper[[i]] <- ifelse(get(paste0("y", i)) > 0, Inf,   0)

      assign(paste0("priors",i), sar_parse_multiple(prior, k, x_obs, i))

      xpx  <- t(x_obs) %*% x_obs       # (X'X)            # k x k
      xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
      xxpxI <- x_obs %*% xpxI          # X(X'X)^(-1)     # n x k (better, compromise)
      AA    <- solve(xpx + solve(get(paste0("priors", i))$diagT)) # (X'X + T^{-1})^{-1}

      # matrices for direct and indirect impacts
      direct[[i]]   <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      indirect[[i]] <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      total[[i]]    <- matrix(NA, ndraw, ncol(x_obs) - 1)    # n x p
      zmean[[i]]    <- matrix(0, n, 1)

      # names of non-constant parameters
      if(cflag == 0) {
        namesNonConstantParams <- colnames(x_obs)
      } else {
        namesNonConstantParams <- colnames(x_obs)[-location.intercept]
      }
      colnames(total[[i]]) <- namesNonConstantParams
      colnames(direct[[i]])   <- namesNonConstantParams
      colnames(indirect[[i]]) <- namesNonConstantParams

      if (computeMarginalEffects) {
        # simulate Monte Carlo estimation of tr(W^i) for i = 1..o before MCMC iterations
        trW.i[[i]] <- tracesWi(W, o=100, iiter=50)
      }

      params <- (k+1) + (J-1)            # parameters beta (k), rho (1), (J-1) cut parameters phi, aber nur (J-2) zu sch?tzen
      # matrix to store the beta + rho parameters for each iteration/draw
      B <- matrix(NA, ndraw, params)
      colnames(B) <- c(paste("beta_", 1:k, sep=""), "rho", paste("y>=", 2:J, sep=""))

      z <- rep(0, n)
      ones <- rep(1, n)
    }
  }

  # truncation points for z, depend only on y, can be precalculated
  #lower <- ifelse(y > 0, 0,  -Inf)
  #upper <- ifelse(y > 0, Inf,   0)
  ####################################################################
  results = c()
  results$nobs = n
  results$nvar = k
  results$x_obs = x_obs
  results$x = x
  results$y = y
  results$W = W
  results$Yk = Yk
  results$cflag = cflag
  results$p = p

  if(nargs() == 5) prior$lflag = 1

  priors = sar_parse(prior,k)

  results$order = priors$order
  results$iter = priors$iter

  output1 = list(); output2 = list(); rmin = list(); rmax = list(); detval = list(); bsave = list()
  for(i in 1:Yk){
    output1[[i]] =  sar_eigs(priors$eflag, get(paste0("W", i)))
    rmin[[i]] = output1[[i]]$rmin
    rmax[[i]] = output1[[i]]$rmax
    output2[[i]] = sar_lndet(priors$ldetflag, get(paste0("W", i)), rmin[[i]], rmax[[i]])
    detval[[i]] = output2[[i]]$detval

    # for every Y, create a matrix for the betas stored in this array
    bsave[[i]] = array(0, dim=c(ndraw*thinning, seq(k,k+Yk-1)[i]))
  }
  rm(i)

  if(priors$mm != 0) rsave = array(rep(0, ndraw*thinning), dim = c(ndraw*thinning, 1, Yk))
  psave = array(rep(0, ndraw*thinning), dim = c(ndraw*thinning, 1, Yk))
  ssave = array(rep(0, ndraw*thinning), dim = c(ndraw*thinning, 1, Yk))
  # vmean = rep(0, n) # this is based on n_obs so this is different so this was defined above in the first Yk loop

  # ====== initializations
  # compute this stuff once to save time
  TI = solve(priors$diagT)
  TIc = TI%*%priors$diffuseprior

  # in.ones = rep(1, n) # this is based on n_obs so this is different so this was defined above in the first Yk loop
  V = in.ones
  vi = in.ones

  # decided to let this be one estimate (corresponding to Y1) instead of Yk estimates
  # Some precalculated quantities for drawing rho from spatialprobit package
  # rho ~ Beta(a1, a2) prior
  lnbprior <- log(beta_prior(detval[[1]][,1], priors$a1, priors$a2))
  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval[[1]]) # do_ldet() gives only 2000 statt 2001 Gridpoints
  # nmk    <- (n-k)/2 # depends on n_obs so defined in the first Yk loop
  rho_grid <- detval[[1]][,1]   # rho grid values
  lndet    <- detval[[1]][,2]   # log-determinant grid values
  rho_gridsq <- rho_grid * rho_grid
  yy       <- (rho_grid[2:nrho] + rho_grid[1:(nrho-1)])

  # define output parameters and lists with priors
  xpy = list(); xpWy = list(); sige = rep(0, Yk); rho = rep(0, Yk); rval = rep(0, Yk);
  ys = list(); Wys = list(); xs = list(); mlike = list(); x_select = list(); TIselect = list(); AA = list()
  sige[seq(1,Yk)] <- priors$sige # defined sige here to make it easier later on
  rho[seq(1,Yk)] <- priors$rho
  rval[seq(1,Yk)] <- priors$rval

  #save(objects,file="fbsTesting.RData")
  #load("fbsTesting.RData")
  #i=1
  #burn.in = 10; ndraw = 100

  if(prior2$novi == 0){ # fit heteroscedastic model
    if(showProgress){pb <- txtProgressBar(min= 1 - burn.in, max=(thinning * ndraw), initial= 1 - burn.in, style=3)}
    iter = 1 - burn.in
    while(iter <= ndraw*thinning){ # start sampling;
      for(i in 1:Yk){
        if(i == 1){
          x_obs = as.matrix(x)   # the imputed Y1 variable is added to the left of the x_obs matrix, so (Yk, Yk-1, ..., Y1, X_obs)
          # x_obs = as.matrix(x_obs[-missing[[i]],]) # remove observations with missing Y_i to get observed cases only
          if(length(missing[[i]]) == 0){ # this should normally not occur
            x_mis = c()
            x_obs = as.matrix(x_obs)
          } else {
            x_mis = t(as.matrix(x_obs[missing[[i]],]))
            #x_obs = as.matrix(x_obs[-missing[[i]],]) # remove observations with missing Y_i to get observed cases only
          }
          xpx = t(x_obs)%*%x_obs
          xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
          xxpxI <- x_obs %*% xpxI              # X(X'X)^(-1)     # n x k (better, compromise)
          Wmm = W[missing[[i]], missing[[i]]] # for the imputation
          Wmo = W[missing[[i]], obs[[i]]] # for the imputation
          assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i))) # this is important because values in y are updated
        } else {
          x_obs = cbind(y[,c(1:seq(0,Yk)[i])], x)   # the imputed Y1 variable is added to the left of the x_obs matrix, so (Yk, Yk-1, ..., Y1, X_obs)
          x_mis = as.matrix(x_obs[missing[[i]],])
          colnames(x_mis) = c(unique(paste0("Y", c(i-1,1))), colnames(x_obs[,-c(1:(i-1))]))
          x_mis = x_mis[,c(i, 1:(i-1), (i+1):ncol(x_mis))]
          #x_obs = as.matrix(x_obs[-missing[[i]],]) # remove observations with missing Y_i to get observed cases only
          colnames(x_obs) = c(unique(paste0("Y", c(i-1,1))), colnames(x_obs[,-c(1:(i-1))]))
          x_obs = as.matrix(x_obs[,c(i, 1:(i-1), (i+1):ncol(x_obs))])
          xpx = t(x_obs)%*%x_obs
          xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
          xxpxI <- x_obs %*% xpxI              # X(X'X)^(-1)     # n x k (better, compromise)
          Wmm = W[missing[[i]], missing[[i]]] # for the imputation
          Wmo = W[missing[[i]], obs[[i]]] # for the imputation
          assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i))) # this is important because values in y are updated
        }

        x_mis = as.matrix(x_mis)
        if(ncol(x_mis) != ncol(x_obs)){x_mis = t(x_mis)} # a small check to make sure x_mis has the right orientation
        assign(paste0("priors",i), sar_parse_multiple(prior, k, x_obs, i))

        # update beta
        xs[[i]] = x_obs * as.double(sqrt(V[[i]]))
        ys[[i]] = sqrt(V[[i]]) * get(paste0("y", i))
        Wys[[i]] = sqrt(V[[i]])* get(paste0("Wy", i))
        k2 = k+(i-1)

        TI = solve(get(paste0("priors",i))$diagT)
        TIc = TI%*%get(paste0("priors",i))$diffuseprior # redefine diffuse prior in sar_parse (priors in array after sarparse?)

        AI = qr.solve(t(xs[[i]])%*%xs[[i]] + sige[i]*TI, diag(k2))
        yss = ys[[i]] - rho[i]*Wys[[i]]
        xpy[[i]] = t(xs[[i]])%*%yss
        b = t(xs[[i]])%*%yss + sige[i]*TIc
        b0 = qr.solve(t(xs[[i]])%*%xs[[i]] + sige[i]*TI, b)
        bhat = norm_rnd(sige[i]*AI) + b0
        xb = xs[[i]]%*%bhat

        # update sige
        n = nrow(xb)
        nu1 = n + 2*get(paste0("priors",i))$nu
        e = (yss - xb)
        d1 = 2*get(paste0("priors",i))$d0 + t(e)%*%e
        chi =  rchisq(1,nu1)
        sige[i] = as.numeric(d1/chi)

        # update vi
        ev = get(paste0("y",i)) - rho[i] * get(paste0("Wy",i)) - x_obs%*%bhat
        chiv = rchisq(n, rval[i]+1)
        # chiv = chi2rnd(rval+1,n,1); % Statistics Toolbox function
        vi[[i]] = ((ev*ev/sige[i]) + in.ones[[i]]*rval[i]) / chiv
        V[[i]] = in.ones[[i]]/vi[[i]]

        # update rval
        if(priors$mm != 0){
          rval[i] = gamm_rnd(1,1) #  get(paste0("priors",i))$mm, get(paste0("priors",i))$kk
        }

        # we use griddy Gibbs to perform rho-draw
        b0 = qr.solve((t(xs[[i]])%*%xs[[i]] + sige[i]*TI), (t(xs[[i]])%*%ys[[i]] + sige[[i]]*TIc))
        bd = qr.solve((t(xs[[i]])%*%xs[[i]] + sige[i]*TI), (t(xs[[i]])%*%Wys[[i]] + sige[[i]]*TIc))
        e0 = ys[[i]] - xs[[i]]%*%b0
        ed = Wys[[i]] - xs[[i]]%*%bd
        epe0 = as.numeric(t(e0)%*%e0)
        eped = as.numeric(t(ed)%*%ed)
        epe0d = as.numeric(t(ed)%*%e0)
        rho[i]  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, rho[i], nmk=nmk[[i]], nrho=nrho, lnbprior, u=u[iter + burn.in])
        if(iter == 1 - burn.in){
          write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        } else {
          write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
        }

        # impute
        #Wy = as.data.frame(get(paste0("Wy", i))[missing[[i]],]) # subset the updated Wy for the missing guys
        #imp.error = rnorm(length(missing[[i]]), 0, sige[i])
        #imp = (rho[[i]] * Wy + x_mis%*%bhat) + imp.error
        #y[missing[[i]], i] <- as.data.frame(imp) # put values in y matrix for imputation of y2 and beyond, and update the estimated values of Wy (above)

        # imputation updated
        Imiss = sparseMatrix(i=1:Nmiss[[i]], j=1:Nmiss[[i]], x=1) # sparse identity matrix
        element1 = qr(Imiss - rho[i]*Wmm)
        element2 = (x_mis%*%bhat) + (rho[i] * Wmo%*%(y[obs[[i]], i]))
        # should be N(0, Imiss*sige[i]) but Imiss is always 1
        imp.error = rnorm(Nmiss[[i]], 0, sige[i])
        element3 = solve(element1, matrix(imp.error))

        # hence:
        imp = as.data.frame(as.matrix(solve(element1, element2) + element3))
        y[missing[[i]], i] = imp[,i]

        if(iter > 0){ # if we are past burn-in, save the draws
          bsave[[i]][iter,] = t(bhat) # bsave[[i]] = Yi, order of coefficients is {Yk, Yk-1, ..., Y1}
          ssave[iter,,i] = sige[i]
          psave[iter,,i] = rho[i]
          vmean[[i]] = vmean[[i]] + vi[[i]]
          y_imp[[i]][iter,] <- imp[,i] # put values in matrix to get posterior conditional distribution of imp values
          if(priors$mm != 0){
            rsave[iter,,i] = rval[i]
          }
        }

        # placed this calculation in the loop at the end because all parameters are available here
        if(iter == ndraw){
          #if(priors$inform_flag == 0){
            xs[[i]] = x_obs * sqrt(V[[i]])
            ys[[i]] = get(paste0("y", i)) * sqrt(V[[i]])
            logdetx = log(det(t(xs[[i]])%*%xs[[i]] + sige[i]*TI))
            mlike[[i]] = sar_marginal_multiple(detval[[i]],e0,ed,epe0,eped,epe0d,nobs=n,k2,logdetx,get(paste0("priors", i))$a1,get(paste0("priors", i))$a2)

          #} else if(priors$inform_flag == 1){
          #  mlike[[i]] = sar_marginal2(detval = detval[[i]],
          #                                      e0 = e0, ed = ed, epe0 = epe0,eped = eped, epe0d = epe0d, nobs = Nobs[[i]], nvar = k2,
          #                                      a1 = priors$a1, a2 = priors$a2, diffuseprior = priors$diffuseprior,TI = TI,xs = xs[[i]],
          #                                      ys = ys[[i]], sige = sige[i], W = get(paste0("W", i)))
          #}
          x_select[[i]] = x_obs
        }
      }
      iter = iter + 1
      if(showProgress){setTxtProgressBar(pb, iter)}
    } # End the sampling
    if(showProgress){close(pb)}
  } else if(prior2$novi == 1){ # fit homoscedastic model
    if(showProgress){pb <- txtProgressBar(min= 1 - burn.in, max=(thinning * ndraw), initial= 1 - burn.in, style=3)}
    iter = 1 - burn.in
    while(iter <= ndraw*thinning){ # start sampling;
      for(i in 1:Yk){
        if(i == 1){
          x_obs = as.matrix(x)   # the imputed Y1 variable is added to the left of the x_obs matrix, so (Yk, Yk-1, ..., Y1, X_obs)
          if(length(missing[[i]]) == 0){ # this should normally not occur
            x_mis = c()
            x_obs = as.matrix(x_obs)
          } else {
            x_mis = t(as.matrix(x_obs[missing[[i]],]))
            #x_obs = as.matrix(x_obs[-missing[[i]],]) # remove observations with missing Y_i to get observed cases only
          }
          xpx = t(x_obs)%*%x_obs
          xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
          xxpxI <- x_obs %*% xpxI              # X(X'X)^(-1)     # n x k (better, compromise)
          Wmm = W[missing[[i]], missing[[i]]] # for the imputation
          Wmo = W[missing[[i]], obs[[i]]] # for the imputation
          assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i))) # this is important because values in y are updated
        } else {
          x_obs = cbind(y[,c(1:seq(0,Yk)[i])], x)   # the imputed Y1 variable is added to the left of the x_obs matrix, so (Yk, Yk-1, ..., Y1, X_obs)
          x_mis = as.matrix(x_obs[missing[[i]],])
          colnames(x_mis) = c(unique(paste0("Y", c(i-1,1))), colnames(x_obs[,-c(1:(i-1))]))
          x_mis = x_mis[,c(i, 1:(i-1), (i+1):ncol(x_mis))]
          #x_obs = as.matrix(x_obs[-missing[[i]],]) # remove observations with missing Y_i to get observed cases only
          colnames(x_obs) = c(unique(paste0("Y", c(i-1,1))), colnames(x_obs[,-c(1:(i-1))]))
          x_obs = as.matrix(x_obs[,c(i, 1:(i-1), (i+1):ncol(x_obs))])
          xpx = t(x_obs)%*%x_obs
          xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
          xxpxI <- x_obs %*% xpxI              # X(X'X)^(-1)     # n x k (better, compromise)
          Wmm = W[missing[[i]], missing[[i]]] # for the imputation
          Wmo = W[missing[[i]], obs[[i]]] # for the imputation
          assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i))) # this is important because values in y are updated
        }

        if(ncol(as.matrix(x_mis)) != ncol(x_obs)){x_mis = t(x_mis)} # a small check to make sure x_mis has the right orientation
        assign(paste0("priors",i), sar_parse_multiple(prior, k, x_obs, i))

        if(method[i] == "probit"){
          if(iter == 1-burn.in){
            TI = solve(get(paste0("priors",i))$diagT)       # Tinv in sarprobit
            TIc = TI%*%get(paste0("priors",i))$diffuseprior # redefine diffuse prior in sar_parse (priors in array after sarparse?)
            TIselect[[i]] = TI
            AA[[i]] = solve(xpx + TIselect[[i]])
          }

          # see LeSage (2009) for choice of burn-in size, often m=5 or m=10 is used!
          # we can also use m=1 together with start.value=z, see LeSage (2009), section 10.1.5
          if (m==1) {
            z = as.double(rtmvnorm.sparseMatrix(n=1, mean=mu[[i]], H=H[[i]], lower=lower[[i]], upper=upper[[i]], burn.in=m, start.value=z))
          } else {
            z = as.double(rtmvnorm.sparseMatrix(n=1, mean=mu[[i]], H=H[[i]], lower=lower[[i]], upper=upper[[i]], burn.in=m))
          }

          # emergency fix if rtmvnorm results in NaN:
          # so replaced that with the normal rtruncnorm function
          #z=c()
          #for(QQ in 1:n){
          #  z[QQ] = rtruncnorm(1, a=lower[[i]][[QQ]], b=upper[[i]][[QQ]], mean = mu[[i]][QQ], sd = 1)
          #}

          # 2. sample from beta | rho, z, y
          Sz <- as.double(S[[i]] %*% z)               # (n x 1); dense
          c1 = (t(x_obs) %*% Sz + TIselect[[i]] %*% get(paste0("priors", i))$diffuseprior)
          c2 <- AA[[i]] %*% c1
          TIselect[[i]] <- AA[[i]]    # no update basically on T, TODO: check this
          beta[[i]] <- as.double(c2 + betadraws[[i]][iter + burn.in, ])

          k2 = dim(xpx)[1]
          bhat = beta[[i]]
          # 3. sample from rho | beta, z
          #---- DRAW RHO ----
          #see LeSage 2009 chapter 5 - page 132 for the explanation of the
          #code below which is used for numerical integration of the rho prior.
          #I changed from the original code to match the notation of the book
          #using c0 and cd below instead of b0 and bd ....
          xpy[[i]] = t(x_obs) %*% z
          xpWy[[i]]= as.double(get(paste0("W", i)) %*% z) # Wz in sarprobit # SW: coerce Wz to vector
          # (from n x 1 sparse matrix! we do not need a sparse matrix here)
          xpWz <- t(x_obs) %*% xpWy[[i]]      # X'Wz      # k x 1
          e0   <- z - xxpxI %*% xpy[[i]]  # z  - X(X'X)^-1X' z
          ed   <- xpWy[[i]] - xxpxI %*% xpWz # Wz - X(X'X)^(-1)X'Wz
          epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
          eped <- as.double(crossprod(ed))
          epe0d<- as.double(crossprod(ed, e0))
          rho[i]  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, rho[i], nmk=nmk[[i]], nrho=nrho, lnbprior, u=u[iter + burn.in])
          if(iter == 1 - burn.in){
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
          } else {
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
          }

          # update S, H and QR decomposition of S and mu after each iteration; before effects
          S[[i]] <- update_I_rW(S[[i]], ind=ind2[[i]], rho[[i]], get(paste0("W", i)))  # update (I - rho * W)
          H[[i]] <- t(S[[i]]) %*% S[[i]]     # H = S'S
          QR <- qr(S[[i]])              # class "sparseQR"

          # impute
          #Wy = as.data.frame(get(paste0("Wy", i))[missing[[i]],]) # subset the updated Wy for the missing guys
          #imp.error = rnorm(length(missing[[i]]), 0, sige[i])
          #imp = (rho[[i]] * Wy + x_mis%*%bhat) + imp.error
          #if(method[i] == "probit") imp = ifelse(imp > mean(mu[[i]]), 1, 0)
          #y[missing[[i]], i] <- as.data.frame(imp) # put values in y matrix for imputation of y2 and beyond, and update the estimated values of Wy (above)

          # imputation updated
          Imiss = sparseMatrix(i=1:Nmiss[[i]], j=1:Nmiss[[i]], x=1) # sparse identity matrix

          element1 = qr(Imiss - rho[i]*Wmm)
          element2 = (x_mis%*%bhat) + (rho[i] * Wmo%*%(y[obs[[i]], i]))
          # should be N(0, Imiss*sige[i]) but Imiss is always 1
          imp.error = rnorm(Nmiss[[i]], 0, sige[i])
          element3 = solve(element1, matrix(imp.error))

          # hence:
          imp = as.data.frame(as.matrix(solve(element1, element2) + element3))
          #if(method[i] == "probit") imp = ifelse(imp > mean(mu[[i]]), 1, 0) # get 0 1 value from z
          y[missing[[i]], i] = imp[,i]

          # update mu
          mu[[i]] <- solve(QR, x_obs %*% bhat)

          # save estimates after burn (i starts negative)
          if (iter > 0) {
            if (thinning == 1) {
              ind[[i]] <- iter
            } else if (iter%%thinning == 0) {
              ind[[i]] <- iter%/%thinning
            } else {
              next
            }

            bsave[[i]][ind[[i]],] = bhat # bsave[[i]] = Yi, order of coefficients is {Yk, Yk-1, ..., Y1}
            zmean[[i]]            = zmean[[i]] + z
            ssave[iter,,i]        = sige[i]
            psave[iter,,i]        = rho[i]
            vmean[[i]]            = vmean[[i]] + vi[[i]]
            y_imp[[i]][iter,]     = imp[,i] # put values in matrix to get posterior conditional distribution of imp values


            # compute effects estimates (direct and indirect impacts) in each MCMC iteration
            if (computeMarginalEffects) {
              o <- 100
              rhovec <- rho[[i]]^(0:(o-1)) # SW: (100 x 1)   mit [1, rho^1, rho^2 ..., rho^99], see LeSage(2009), eqn (4.145), p.115
              if( cflag == 1 ){ #has intercept
                location.intercept = which(colnames(x_obs) == "Intercept")
                beff <- beta[[i]][-location.intercept] # beff is parameter vector without constant
              }else if(cflag == 0){
                beff <- beta[[i]]     # no constant in model
              }
              # beff is parameter vector without constant!
              # See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
              #   direct: M_r(D) = n^{-1} tr(S_r(W))           # SW: efficient approaches available, see chapter 4, pp.114/115
              #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # SW: Problem: S_r(W) is dense, but can be solved via QR decomposition of S
              # indirect: M_r(I) = M_r(T) - M_r(D)
              # SW: See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
              pdfz <- dnorm(as.numeric(mu[[i]]))                     # standard normal pdf phi(mu)
              dd   <- sparseMatrix(i=1:n, j=1:n, x=pdfz)       # dd is diagonal matrix with pdfz as diagonal (n x n)

              dir      <- as.double(t(pdfz) %*% trW.i[[i]] %*% rhovec / n)  # (1 x n) * (n x o) * (o x 1)
              # direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * beta_r
              avg_direct     <- dir * beff      # (p x 1)

              # We compute the average total effects without inverting S
              # unlike in the LeSage Matlab Code,
              # but using the QR decomposition of S which we already have!
              # average total effects = n^(-1) * 1_n' %*% (D %*% S^(-1) * b[r]) %*% 1_n
              #                       = n^(-1) * 1_n' %*% (D %*% x) * b[r]
              # where D=dd is the diagonal matrix containing phi(mu)
              # and x is the solution of S %*% x = 1_n, obtained from the QR-decompositon
              # of S. The average total effects is then the mean of (D %*% x) * b[r]
              # average total effects, which can be furthermore done for all b[r] in one operation.
              avg_total    <- mean(dd %*% qr.coef(QR, in.ones[[i]])) * beff
              avg_indirect <- avg_total - avg_direct    # (p x 1)

              total[[i]][ind[[i]], ]      <- avg_total    # an (ndraw-nomit x p) matrix
              direct[[i]][ind[[i]], ]     <- avg_direct   # an (ndraw-nomit x p) matrix
              indirect[[i]][ind[[i]], ]   <- avg_indirect # an (ndraw-nomit x p) matrix
            }
            x_select[[i]] = x_obs
          }

        } else if(method[i] == "orderedprobit"){

          if(iter == 1-burn.in){
            TI = solve(get(paste0("priors",i))$diagT)       # Tinv in sarprobit
            TIc = TI%*%get(paste0("priors",i))$diffuseprior # redefine diffuse prior in sar_parse (priors in array after sarparse?)
            TIselect[[i]] = TI
            AA[[i]] = solve(xpx + TIselect[[i]])
          }

          # determine lower and upper bounds for z depending on the value of y and phi
          # eqn (10.14), p.299
          # besser: y = 1..J
          lower <- phi[get(paste0("y", i))]
          upper <- phi[get(paste0("y", i)) + 1]

          # see LeSage (2009) for choice of burn-in size, often m=5 or m=10 is used!
          # we can also use m=1 together with start.value=z, see LeSage (2009), section 10.1.5
          if (m==1) {
            z = as.double(rtmvnorm.sparseMatrix(n=1, mean=mu[[i]], H=H[[i]], lower=lower[[i]], upper=upper[[i]], burn.in=m, start.value=z))
          } else {
            z = as.double(rtmvnorm.sparseMatrix(n=1, mean=mu[[i]], H=H[[i]], lower=lower[[i]], upper=upper[[i]], burn.in=m))
          }

          # emergency fix if rtmvnorm results in NaN:
          # so replaced that with the normal rtruncnorm function
          #z=c()
          #for(QQ in 1:n){
          #  z[QQ] = rtruncnorm(1, a=lower[[i]][[QQ]], b=upper[[i]][[QQ]], mean = mu[[i]][QQ], sd = 1)
          #}
          # test
          #test = as.data.frame(cbind(z, y[-5,]))
          #plot(test[test$V2 == 1, ]$z) # up to 6 for year of study

          # 2. sample from beta | rho, z, y
          c_prior <- AA[[i]]  %*% (t(x_obs) %*% S[[i]] %*% z + TI %*% c_prior)
          TIselect[[i]] <- AA[[i]] # no update basically on T, TODO: check this
          beta[[i]] <- as.vector(rmvnorm(n=1, mean=c_prior, sigma=TIselect[[i]]))

          k2 = dim(xpx)[1]
          bhat = beta[[i]]
          # 3. sample from rho | beta, z
          #---- DRAW RHO ----
          #see LeSage 2009 chapter 5 - page 132 for the explanation of the
          #code below which is used for numerical integration of the rho prior.
          #I changed from the original code to match the notation of the book
          #using c0 and cd below instead of b0 and bd ....
          xpy[[i]] = t(x_obs) %*% z
          xpWy[[i]]= as.double(get(paste0("W", i)) %*% z) # Wz in sarprobit # SW: coerce Wz to vector
          # (from n x 1 sparse matrix! we do not need a sparse matrix here)
          xpWz <- t(x_obs) %*% xpWy[[i]]      # X'Wz      # k x 1
          e0   <- z - xxpxI %*% xpy[[i]]  # z  - X(X'X)^-1X' z
          ed   <- xpWy[[i]] - xxpxI %*% xpWz # Wz - X(X'X)^(-1)X'Wz
          epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
          eped <- as.double(crossprod(ed))
          epe0d<- as.double(crossprod(ed, e0))
          rho[i]  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, rho[i], nmk=nmk[[i]], nrho=nrho, lnbprior, u=u[iter + burn.in])
          if(iter == 1 - burn.in){
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
          } else {
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
          }

          # 4. determine bounds/cut-points p(phi_j | phi_{-j}, z, y, beta) for j = 2,...,J-1
          # phi_j = 0 is set fixed!
          for (j in 2:(J-1)) {
            phi.lower <- max(max(z[get(paste0("y", i)) == j]),        phi[j-1+1])   # \bar{phi}_{j-1}, SW: +1 is needed as our vector index starts with 1
            phi.upper <- min(min(z[get(paste0("y", i)) == (j + 1) ]), phi[j+1+1])   # \bar{phi}_{j+1}

            # Sample phi_{j | phi_{-j}, z, y, beta)
            phi[j + 1]   <- runif(n=1, min=phi.lower, max=phi.upper)
          }

          # update S, H and QR decomposition of S and mu after each iteration; before effects
          S[[i]] <- update_I_rW(S=S[[i]], ind=ind2[[i]], rho=rho[i], W=get(paste0("W", i)))  # update (I - rho * W)
          H[[i]] <- t(S[[i]]) %*% S[[i]]     # H = S'S
          QR <- qr(S[[i]])              # class "sparseQR"

          # impute
          #Wy = as.data.frame(get(paste0("Wy", i))[missing[[i]],]) # subset the updated Wy for the missing guys
          #imp.error = rnorm(length(missing[[i]]), 0, sige[i])
          #imp = (rho[[i]] * Wy + x_mis%*%bhat) + imp.error
          #if(method[i] == "orderedprobit"){
          #  imp <- cut(as.double(imp), breaks=phi, labels=FALSE, ordered_result = TRUE) # split according to phi values
          #}
          #y[missing[[i]], i] <- as.data.frame(imp) # put values in y matrix for imputation of y2 and beyond, and update the estimated values of Wy (above)

          # imputation updated
          Imiss = sparseMatrix(i=1:Nmiss[[i]], j=1:Nmiss[[i]], x=1) # sparse identity matrix
          element1 = qr(Imiss - rho[i]*Wmm)
          element2 = (x_mis%*%bhat) + (rho[i] * Wmo%*%(y[obs[[i]], i]))
          # should be N(0, Imiss*sige[i]) but Imiss is always 1
          imp.error = rnorm(Nmiss[[i]], 0, sige[i])
          element3 = solve(element1, matrix(imp.error))

          # hence:
          imp = as.data.frame(as.matrix(solve(element1, element2) + element3))
          if(method[i] == "orderedprobit"){
            imp <- cut(as.double(imp), breaks=phi, labels=FALSE, ordered_result = TRUE) # split according to phi values
          }
          y[missing[[i]], i] = imp[,i]

          # update mu
          mu[[i]] <- solve(QR, x_obs %*% bhat)

          # save estimates after burn (i starts negative)
          if (iter > 0) {
            if (thinning == 1) {
              ind[[i]] <- iter
            } else if (iter%%thinning == 0) {
              ind[[i]] <- iter%/%thinning
            } else {
              next
            }

            bsave[[i]][ind[[i]],] = bhat # bsave[[i]] = Yi, order of coefficients is {Yk, Yk-1, ..., Y1}
            zmean[[i]]            = zmean[[i]] + z
            ssave[iter,,i]        = sige[i]
            psave[iter,,i]        = rho[i]
            vmean[[i]]            = vmean[[i]] + vi[[i]]
            y_imp[[i]][iter,]     = imp[,i] # put values in matrix to get posterior conditional distribution of imp values
            B[iter,] <- c(beta[[i]], rho[[i]], phi[2:J])   # (k + 1) + (J - 1)

            x_select[[i]] = x_obs
          }
        } else {

          xpy[[i]] = t(x_obs) %*% get(paste0("y", i))
          xpWy[[i]] = t(x_obs) %*% get(paste0("Wy", i))

          assign(paste0("priors",i), sar_parse_multiple(prior, k, x_obs, i))
          TI = solve(get(paste0("priors",i))$diagT)
          TIc = TI %*% get(paste0("priors",i))$diffuseprior # redefine diffuse prior in sar_parse (priors in array after sarparse?)

          k2 = dim(xpx)[1]

          AI = qr.solve(xpx + sige[i]*TI, diag(k2))
          ys = get(paste0("y", i)) - rho[i]*get(paste0("Wy", i))
          b = t(as.matrix(x_obs)) %*% ys + sige[i]*TIc
          b0 = qr.solve(xpx + sige[i]*TI, b)
          bhat = norm_rnd(sige[i]*AI) + b0
          xb = x_obs %*% bhat

          # update sige
          n = nrow(xb)
          nu1 = n + 2*get(paste0("priors",i))$nu
          e = (ys - xb)
          d1 = 2*get(paste0("priors",i))$d0 + t(e)%*%e
          chi =  rchisq(1,nu1)
          sige[i] = as.numeric(d1/chi)

          # update rho using griddy Gibbs
          AI = qr.solve(xpx + sige[i]*TI, diag(k2))
          b0 = qr.solve(xpx + sige[i]*TI, xpy[[i]] + sige[i]*TIc)
          bd = qr.solve(xpx + sige[i]*TI, xpWy[[i]] + sige[i]*TIc)
          e0 = get(paste0("y",i)) - as.matrix(x_obs)%*%b0
          ed = get(paste0("Wy",i)) - as.matrix(x_obs)%*%bd
          epe0 = as.numeric(t(e0) %*% e0)
          eped = as.numeric(t(ed) %*% ed)
          epe0d = as.numeric(t(ed) %*% e0)
          rho[i]  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, rho[i], nmk=nmk[[i]], nrho=nrho, lnbprior, u=u[iter + burn.in])
          if(iter == 1 - burn.in){
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
          } else {
            write.table(rho[i], paste0("rhoestimate", model, "fbs.txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
          }

          # impute
          #Wy = as.data.frame(get(paste0("Wy", i))[missing[[i]],]) # subset the updated Wy for the missing guys
          #imp.error = rnorm(length(missing[[i]]), 0, sige[i])
          #imp = (rho[[i]] * Wy + x_mis%*%bhat) + imp.error
          #y[missing[[i]], i] <- as.data.frame(imp) # put values in y matrix for imputation of y2 and beyond

          # imputation updated
          Imiss = sparseMatrix(i=1:Nmiss[[i]], j=1:Nmiss[[i]], x=1) # sparse identity matrix

          element1 = qr(Imiss - rho[i]*Wmm)
          element2 = (x_mis%*%bhat) + (rho[i] * Wmo%*%(y[obs[[i]], i]))
          # should be N(0, Imiss*sige[i]) but Imiss is always 1
          imp.error = rnorm(Nmiss[[i]], 0, sige[i])
          element3 = solve(element1, matrix(imp.error))

          # hence:
          imp = as.data.frame(as.matrix(solve(element1, element2) + element3))
          y[missing[[i]], i] = imp[,i]

          if(iter > 0){ # if we are past burn-in, save the draws
            bsave[[i]][iter,] = t(bhat) # bsave[[i]] = Yi, order of coefficients is {Yk, Yk-1, ..., Y1}
            ssave[iter,,i] = sige[i]
            psave[iter,,i] = rho[i]
            vmean[[i]] = vmean[[i]] + vi[[i]]
            y_imp[[i]][iter,] <- imp[,i] # put values in matrix to get posterior conditional distribution of imp values
          }

          x_select[[i]] = x_obs
        }
      }
      iter = iter + 1
      if(showProgress){setTxtProgressBar(pb, iter)}
    } # End the sampling
    if(showProgress){close(pb)}

  } else {
    stop('unrecognized prior2.novi value on input; cannot decide whether to fit a homo- or heteroscedastic model')
  }

  # pre-calculate traces for the x-impacts calculations
  uiter = 50
  maxorderu = 100
  rv = list()
  for(i in 1:Yk){
    rv[[i]] = matrix(rnorm(n * uiter), n, uiter)
  }
  tracew = matrix(0, maxorderu, Yk)
  wjjju = rv
  for(i in 1:Yk){
    for(jjj in 1:maxorderu){
      wjjju[[i]] = get(paste0("W", i))%*%wjjju[[i]]
      tracew[jjj,i] = mean(mean(rv[[i]]*wjjju[[i]]))
    }
  }

  traces = tracew
  traces[1,] = 0
  bdraws = list()
  for(i in 1:Yk){
    if(method[i] != "probit"){
      traces[2,i] = sum(sum(t(get(paste0("W", i)))*get(paste0("W", i))))/n
      if(cflag == 1){
        bdraws[[i]] = bsave[[i]][,c(2:ncol(bsave[[i]]))]
      } else if(cflag == 0){
        bdraws = bsave
      }

      total[[i]] = array(0, dim=c(ndraw*thinning, seq(p,p+Yk-1)[i]))
      direct[[i]] = array(0, dim=c(ndraw*thinning, seq(p,p+Yk-1)[i]))
      indirect[[i]] = array(0, dim=c(ndraw*thinning, seq(p,p+Yk-1)[i]))
    }
  }
  trs = rbind(1, traces)
  ntrs = nrow(trs)
  trbig = t(trs)

  pdraws = psave

  ree = seq(0, ntrs-1)
  rmat = rep(0, ntrs)

  #for(z in 1:Yk){
  #  for(i in 1:(ndraw-burn.in)){
  #    rmat = pdraws[i,,z]^ree
  #    for(j in 1:p){
  #      beta = bdraws[[z]][i,j]
  #      brmat = beta*rmat
  #      btrmat = (beta*trbig)*rmat
  #      for(k in 1:ntrs){        <- has to be 3 of the array, which is a problem (4 dimensional array needed)
  #        total[i,j,k] = brmat[k]
  #        direct[i,j,k] = btrmat[k]
  #        indirect[i,j,k] = total[i,j,k] - direct[i,j,k]
  #      }
  #    }
  #  }
  #}

  # compute posterior means and log marginal likelihood for return arguments
  beta = list()
  yhat = list()
  e = list()
  epe = rep(0,Yk)
  ym = list()
  rsqr2 = rep(0, Yk)
  beta_std = list(); sige_std = rep(0, Yk); rho_std = rep(0, Yk)
  for(i in 1:Yk){
    beta[[i]] = colMeans(bsave[[i]])
    rho[i] = mean(psave[,,i])
    sige[i] = mean(ssave[,,i])
    vmean[[i]] = vmean[[i]]/(ndraw)
    V[[i]] = in.ones[[i]]/vmean[[i]]

    # compute R-squared
    yhat[[i]] = qr.solve(Diagonal(n) - rho[i]*get(paste0("W", i)), x_select[[i]]%*%beta[[i]])
    e[[i]] = get(paste0("y", i)) - yhat[[i]]
    epe[i] = t(e[[i]])%*%e[[i]]
    ym[[i]] = get(paste0("y", i)) - mean(get(paste0("y", i)))
    rsqr2[i] = t(ym[[i]])%*%ym[[i]]
    beta_std[[i]] = apply(bsave[[i]], 2, sd)
    sige_std[i] = sd(ssave[,,i])
    rho_std[i] = sd(psave[,,i])
  }

  if("orderedprobit" %in% method){
    phi   <- c(-Inf, colMeans(B)[(k+2):((k+2) + (J-2))], Inf)
    fitted.values   <- solve(qr(S[[1]]), x_obs %*% beta[[1]]) # z = (I_n - rho * W)^{-1}(X * beta)
    fitted.response <- cut(as.double(fitted.values), breaks=phi, labels=FALSE, ordered_result = TRUE) # split according to phi values
  }

  results$sige = sige

  # compute R-squared
  sige = epe / as.numeric(n) - seq(k, k + (i-1))
  results$sigma = sige
  rsqr1 = epe
  results$rsqr = 1 - rsqr1 / rsqr2     # r-squared
  rsqr1 = rsqr1 / as.numeric(n) - seq(k, k + (i-1))
  rsqr2 = rsqr2 / (as.numeric(n) - 1)
  results$rbar = 1 - (rsqr1/rsqr2)  # rbar-squared

  # write output
  results$meth  = 'sar_g'
  results$total = total
  results$direct = direct
  results$indirect = indirect
  results$beta_std = beta_std
  results$sige_std = sige_std
  results$rho_std = rho_std
  results$beta = beta
  results$rho = rho
  results$y_imp = y_imp
  results$bdraw = bsave
  results$pdraw = psave
  results$sdraw = ssave
  results$mlike = mlike
  results$vmean = vmean
  results$zmean = zmean
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



