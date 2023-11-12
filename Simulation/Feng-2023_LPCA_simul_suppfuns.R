#############################################
##### Large Dim. Nonlinear Factor Model #####
### Define supporting functions #############
######## for simulation #####################
######## date: Nov 12, 2023 #################
#############################################
# large-dim x
hdfun1 <- function(u, v, sig=.1) {
  out <- (1/sqrt(2*pi)/sig) * exp(-(u-v)^2/sig^2)
  return(out)
}

hdfun2 <- function(u, v, sig=.1) {
  out <- exp(-abs(u-v)/sig)
  return(out)
}

hdfun3 <- function(u,v) {
  out <- 1-(1+exp(15*(0.8*abs(u-v))^0.8-0.1))^(-1)
  return(out)
}

funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2, hdfun3=hdfun3)


# DGP generator
dgp <- function(n, p, hdmodel, sigma=1, binary=F) {
    # latent variable
    alpha <- runif(n, 0, 1)
    eval.ind <- sapply(round(c(0.1, 0.5, 0.9)*n), function(i) nth(alpha, i, index.return=T))
    
    # HD covariates
    eta  <- runif(p, 0, 1)
    mean <- outer(alpha, eta, FUN = hdmodel) 
  
    if (binary) {
      x    <- 1 * (matrix(runif(n*p,0,1), n, p) <= mean)
    } else {
      x    <- mean + matrix(rnorm(n*p, 0, sigma), n, p)   # n by p matrix
    }
    
    return(list(x=x, mean=mean, eval.ind=eval.ind))
}

################################################
# Note the following KNN can accept a vector of Ks and report a list of matrices
# Searching for KNN
findknn <- function(v, A2, Kseq=NULL) {
     tmp   <- abs(A2 - v[-1])
     diag(tmp)  <- -1        # implicitly delete diagonal elements
     tmp[v[1],] <- -1        # implicitly delete the v[1]th element
     tmp.d <- colMaxs(tmp, value = T)
     out   <- sapply(Kseq, function(K) ifelse(tmp.d <= nth(tmp.d, K), TRUE, FALSE)) # n by length(Kseq) matrix
     return(out)
}

knn.index <- function(A, Kseq=NULL) {    # A: n by p; # Kseq: vector of tuning param.
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p
    Kmat <- apply(cbind(1:nrow(A), A2), 1, function(v) findknn(v=v, A2=A2, Kseq=Kseq))
    return(Kmat)           # a matrix: neighborhood for each unit saved in each column, each n rows -> one choice of K in Kseq 
}


# Local PCA (for each point)
# A: n by p, input matrix; index: n by 1, KNN index; nlam: MAXimum number of vectors to be extracted; i: unit of interest
lpca <- function(A, index, nlam, n, K, i) {      
    A <- cbind(1:n, A)   # add a row number
    A <- A[index,]
    svd <- irlba(A[,-1], nu=nlam, nv=nlam)
    pos <- which(A[,1]==i)
    r <- sel.r(svd$d, K)
    
    # multiply by singular values
    if (r==1) {
      mean  <- tcrossprod(svd$u[,1], (svd$v[,1]*svd$d[1]))[pos,]
    } else {
      mean <- tcrossprod(svd$u[,1:r], (svd$v[,1:r] %*% diag(svd$d[1:r])))[pos,]
    }

    return(c(mean))    # a vector of length p*2  [lpca.nlam]  
}

# Prediction (output: mean prediction for each i with a seq of K's)
pred <- function(data, fullkmat, nlam, n, Kseq, i) {
  L <- nrow(fullkmat)/n
  out <- sapply(c(1:L), function(l) lpca(A=data, index=fullkmat[((l-1)*n+1):(l*n),i], nlam=nlam, n=n, K=Kseq[l], i=i))   # (p/2 * 1) by length(Kseq)
  
  return(out)
}

# Computation, for all units in range
compute <- function(range, data, Kseq, nlam, n, p, prop.mat=0.5) {
    p1 <- round(p*prop.mat)
    fullkmat <- knn.index(A=data[,1:p1], Kseq=Kseq)
    x <- data[,(p1+1):p]
    fit  <- sapply(range, function(i) pred(data=x, fullkmat=fullkmat, nlam=nlam, n=n, Kseq=Kseq, i=i))      
    
    return(t(fit))     # length(range) by [(p/2 *1) * length(kseq)]     (lpca.nlam) * Kseq
}

# compute statistics
comstat <- function(est, true, ind) {
  dev <- abs(est - true)
  return(c(max(dev), mean(dev), dev[ind, ncol(est)]))
}


# sim function
sim <- function(i, n, p, model, Kseq, nlam, sigma=1, binary=F, prop.mat=0.5) {
  p1 <- round(p*prop.mat); p2 <- p-p1
  # generate data
  data <- dgp(n=n, p=p, hdmodel=model, sigma=sigma, binary=binary)
  mean <- data$mean[,(p1+1):p]
  x <- data$x; x[data$eval.ind, p] <- 0          # make three units in the last col. = 0 
  
  # LPCA
  range <- 1:n
  fit   <- compute(range, data=x, Kseq, nlam, n, p, prop.mat)
  
  # return 5 by length(Kseq)
  g <- 1
  out <- sapply(c(1:(g*length(Kseq))), function(l) comstat(fit[,((l-1)*p2+1):(l*p2)], mean, data$eval.ind))
  
  # Linear PCA
  x.gpca <- x
  # select r based on demeaned data
  col.mean <- colmeans(x.gpca); x.gpca1 <- sweep(x.gpca, 2, col.mean)
  est.r <- PCA_FN(x.gpca1, 9)
  fit.linear  <- irlba(x.gpca1, nu=est.r, nv=est.r)
  if (est.r==1) {
    mean.linear <- sweep(tcrossprod(fit.linear$u, fit.linear$v * fit.linear$d), 2, col.mean, "+")
  } else {
    mean.linear <- sweep(tcrossprod(fit.linear$u, fit.linear$v %*% diag(fit.linear$d)), 2, col.mean, "+") 
  }
  stat.linear <- comstat(mean.linear[,(p1+1):p], mean, data$eval.ind)
  out         <- cbind(out, stat.linear, rep(est.r+1, length(stat.linear)))

  # double demeaned data
  row.mean <- rowmeans(x.gpca1); x.gpca2 <- x.gpca1 - row.mean
  est.r <- PCA_FN(x.gpca2, 9)
  fit.linear  <- irlba(x.gpca2, nu=est.r, nv=est.r)
  if (est.r==1) {
    mean.linear <- sweep(tcrossprod(fit.linear$u, fit.linear$v * fit.linear$d) + row.mean, 2, col.mean, "+")
  } else {
    mean.linear <- sweep(tcrossprod(fit.linear$u, fit.linear$v %*% diag(fit.linear$d)) + row.mean, 2, col.mean, "+") 
  }
  stat.linear <- comstat(mean.linear[,(p1+1):p], mean, data$eval.ind)
  out         <- cbind(out, stat.linear, rep(est.r+2, length(stat.linear)))
  
  return(out)   # 5 by (Kseq * 1 + 2 + 2)  [(lpca.nlam) * Kseq, linear, linear.demean]
}

##################################################
sel.r <- function(sval, n) {
  d <- length(sval)
  ratio <- (sval[-d] / sval[-1]) < log(log(n))
  if (any(ratio)) {
    r <- max(which.max(ratio)-1, 1)
  } else {
    r <- d - 1
  }
  return(r)
}

