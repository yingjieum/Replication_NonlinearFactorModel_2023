#############################################
##### Large Dim. Nonlinear Factor Model #####
### Define supporting functions #############
######## for simulation #####################
######## date: Nov 12, 2023 #################
#############################################

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

sc.pred <- function(data, ind, post, eq=1, lb=0) {
  Z <- t(data[-ind, -post])
  d.n <- ncol(Z); d.p <- nrow(Z)
  
  fit <- NULL
  for (i in ind) {
    x <- Variable(d.n)
    
    # define a quad objective fn
    obj <- Minimize(quad_form(data[i,-post] - Z %*% x, diag(d.p)))
    constraints <- list(sum_entries(x) == eq, x >= lb)
    prob <- Problem(obj, constraints)
    sol  <- solve(prob, solver="ECOS")
    
    w      <- sol$getValue(x)
    fit <- cbind(fit, t(data[-ind,]) %*% w)
  }
  
  return(fit)
}
