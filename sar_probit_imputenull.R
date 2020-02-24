sar_probit_mcmc = function (y, X, W, ndraw = 1000, burn.in = 100, thinning = 1, 
          prior = list(a1 = 1, a2 = 1, c = rep(0, ncol(X)), T = diag(ncol(X)) * 1e+12, lflag = 0),
          start = list(rho = 0.75, beta = rep(0, ncol(X))), m = 10, 
          computeMarginalEffects = TRUE, showProgress = FALSE){
  
  #X = as.matrix(X)
  #ndraw = 1000
  #burn.in = 100
  #thinning = 1
  #prior = list(a1 = 1, a2 = 1, c = rep(0, ncol(X)), T = diag(ncol(X)) * 1e+12, lflag = 0)
  #start = list(rho = 0.75, beta = rep(0, ncol(X)))
  #m = 10
  #computeMarginalEffects = TRUE
  #showProgress = FALSE
  
  timet <- Sys.time()
  n <- nrow(X)
  n1 <- nrow(X)
  n2 <- nrow(W)
  k <- ncol(X)
  I_n <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
  if (is.null(colnames(X))) 
    colnames(X) <- paste("x", 1:k, sep = "")
  if (length(c(which(y == 0), which(y == 1))) != length(y)) {
    stop("sarprobit: not all y-values are 0 or 1")
  }
  if (n1 != n2 && n1 != n) {
    stop("sarprobit: wrong size of spatial weight matrix W")
  }
  if (!inherits(W, "sparseMatrix") || any(diag(W) != 0)) {
    stop("sarprobit: spatial weights matrix W must be a sparse matrix with zeros in the main diagonal")
  }
  ind <- match(n, apply(X, 2, sum))
  if (is.na(ind)) {
    cflag <- 0
    p <- k
  } else if (ind == 1) {
    cflag <- 1
    p <- k - 1
  } else {
    stop("sarprobit: intercept term must be in first column of the X-matrix")
  }
  rho <- start$rho
  beta <- start$beta
  c <- rep(0, k)
  if (is.numeric(prior$c) && length(prior$c) == k) {
    c <- prior$c
  }
  if (is.matrix(prior$T) && ncol(prior$T) == k && isSymmetric(prior$T) && 
      det(prior$T) > 0) {
    T <- prior$T
  } else {
    T <- diag(k) * 1e+12
  }
  Tinv <- solve(T)
  if (class(W) == "dgCMatrix") {
    I <- sparseMatrix(i = 1:n, j = 1:n, x = Inf)
    S <- (I - rho * W)
    ind <- which(is.infinite(S@x))
    ind2 <- which(!is.infinite(S@x))
    S@x[ind] <- 1
  } else {
    S <- I_n - rho * W
  }
  H <- t(S) %*% S
  QR <- qr(S)
  mu <- solve(QR, X %*% beta)
  lower <- ifelse(y > 0, 0, -Inf)
  upper <- ifelse(y > 0, Inf, 0)
  rmin <- -1
  rmax <- 1
  lflag <- 0
  if (is.numeric(prior$lflag) && lflag %in% c(0, 1, 2)) 
    lflag <- prior$lflag
  tmp <- sar_lndet(lflag, W, rmin, rmax)
  detval <- tmp$detval
  a1 <- 1
  a2 <- 1
  if (is.numeric(prior$a1)) 
    a1 <- prior$a1
  if (is.numeric(prior$a2)) 
    a2 <- prior$a2
  lnbprior <- log(beta_prior(detval[, 1], a1, a2))
  u <- runif(thinning * ndraw + burn.in)
  nrho <- nrow(detval)
  nmk <- (n - k)/2
  rho_grid <- detval[, 1]
  lndet <- detval[, 2]
  rho_gridsq <- rho_grid * rho_grid
  yy <- (rho_grid[2:nrho] + rho_grid[1:(nrho - 1)])
  B <- matrix(NA, ndraw, k + 1)
  if (showProgress) {
    pb <- txtProgressBar(min = 0, max = (thinning * ndraw + 
                                           burn.in), initial = 0, style = 3)
  }
  tX <- t(X)
  xpx <- t(X) %*% X
  xpxI <- solve(xpx)
  xxpxI <- X %*% xpxI
  AA <- solve(xpx + Tinv)
  betadraws <- rmvnorm(n = (burn.in + ndraw * thinning), mean = rep(0, 
                                                                    k), sigma = AA)
  direct <- matrix(NA, ndraw, p)
  indirect <- matrix(NA, ndraw, p)
  total <- matrix(NA, ndraw, p)
  zmean <- rep(0, n)
  if (cflag == 0) {
    namesNonConstantParams <- colnames(X)
  }  else {
    namesNonConstantParams <- colnames(X)[-1]
  }
  colnames(total) <- namesNonConstantParams
  colnames(direct) <- namesNonConstantParams
  colnames(indirect) <- namesNonConstantParams
  if (computeMarginalEffects) {
    trW.i <- tracesWi(W, o = 100, iiter = 50)
  }
  z <- rep(0, n)
  ones <- rep(1, n)
  
  y_imp = matrix(0, ndraw*thinning, n)
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
    if (m == 1) {
      #z <- as.double(rtmvnorm.sparseMatrix(n = 1, mean = mu, 
      #                                     H = H, lower = lower, upper = upper, burn.in = m, 
      #                                     start.value = z))
      z=c()          
      for(QQ in 1:n){
        z[QQ] = rtruncnorm(1, a=lower[[QQ]], b=upper[[QQ]], mean = mu[QQ], sd = 1)
      }
    } else {
      #z <- as.double(rtmvnorm.sparseMatrix(n = 1, mean = mu, 
      #                                     H = H, lower = lower, upper = upper, burn.in = m))
      z=c()          
      for(QQ in 1:n){z[QQ] = rtruncnorm(1, a=lower[[QQ]], b=upper[[QQ]], mean = mu[QQ], sd = 1)}
    }
    Sz <- as.double(S %*% z)
    c2 <- AA %*% (tX %*% Sz + Tinv %*% c)
    T <- AA
    beta <- as.double(c2 + betadraws[i + burn.in, ])
    xpz <- tX %*% z
    Wz <- as.double(W %*% z)
    xpWz <- tX %*% Wz
    e0 <- z - xxpxI %*% xpz
    ed <- Wz - xxpxI %*% xpWz
    epe0 <- as.double(crossprod(e0))
    eped <- as.double(crossprod(ed))
    epe0d <- as.double(crossprod(ed, e0))
    rho <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, 
                    eped, epe0d, rho, nmk = nmk, nrho = nrho, lnbprior, 
                    u = u[i + burn.in])
    S <- update_I_rW(S, ind = ind2, rho, W)
    H <- t(S) %*% S
    QR <- qr(S)
    mu <- solve(QR, X %*% beta)
    
    write.table(rho, paste0("rhoestimate_nullmodel_lflag", lflag, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = TRUE)

    # save the output
    Imiss = I_n # sparse identity matrix
    element1 = qr(Imiss - rho*W)
    element2 = (X%*%beta) + (rho * W%*%y)
    imp.error = rnorm(n, 0, 1)
    element3 = solve(element1, matrix(imp.error))
    
    # hence:
    imp = solve(element1, element2) + element3

    if (i > 0) {
      y_imp[i,] = t(as.matrix(imp)) # put values in matrix to get posterior conditional distribution of imp values
      if (thinning == 1) {
        ind <- i
      }
      else if (i%%thinning == 0) {
        ind <- i%/%thinning
      }
      else {
        next
      }
      B[ind, ] <- c(beta, rho)
      zmean <- zmean + z
      if (computeMarginalEffects) {
        o <- 100
        rhovec <- rho^(0:(o - 1))
        if (cflag == 1) {
          beff <- beta[-1]
        }
        else if (cflag == 0) {
          beff <- beta
        }
        pdfz <- dnorm(as.numeric(mu))
        dd <- sparseMatrix(i = 1:n, j = 1:n, x = pdfz)
        dir <- as.double(t(pdfz) %*% trW.i %*% rhovec/n)
        avg_direct <- dir * beff
        avg_total <- mean(dd %*% qr.coef(QR, ones)) * 
          beff
        avg_indirect <- avg_total - avg_direct
        total[ind, ] <- avg_total
        direct[ind, ] <- avg_direct
        indirect[ind, ] <- avg_indirect
      }
    }
    if (showProgress) 
      setTxtProgressBar(pb, i + burn.in)
  }
  if (showProgress) 
    close(pb)
  beta <- colMeans(B)[1:k]
  rho <- colMeans(B)[k + 1]
  S <- (I_n - rho * W)
  fitted.values <- solve(qr(S), X %*% beta)
  fitted.response <- as.numeric(fitted.values >= 0)
  results <- NULL
  results$time <- Sys.time() - timet
  results$nobs <- n
  results$nvar <- k
  results$y <- y
  results$zip <- n - sum(y)
  results$beta <- colMeans(B)[1:k]
  results$rho <- colMeans(B)[k + 1]
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values
  results$fitted.response <- fitted.response
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1 <- a1
  results$a2 <- a2
  results$rmax <- rmax
  results$rmin <- rmin
  results$tflag <- "plevel"
  results$lflag <- lflag
  results$cflag <- cflag
  results$lndet <- detval
  results$names <- c(colnames(X), "rho")
  results$B <- B
  results$bdraw <- B[, 1:k]
  results$pdraw <- B[, k + 1]
  results$total <- total
  results$direct <- direct
  results$indirect <- indirect
  results$imputed <- y_imp
  results$W <- W
  results$X <- X
  class(results) <- "sarprobit"
  return(results)
}
