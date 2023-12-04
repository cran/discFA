# -log likelihood for a Poisson sample.
dpoliku <- function(x){
  t <- mean(x)    # The MLE equals the sample mean.
  f <- -sum(log(exp(-t) * t^x / factorial(x)))
  return(f)
}

# -log likelihood for dependent Poisson in dimension d.
dpolikgq <- function(t0,x){
  x <- as.matrix(x)
  t <- colMeans(x) - t0 # This is where we use the parameter restriction."
  n <- nrow(x)
  d <- ncol(x)
  liksum <- 0
  # The following loop calculates the log likelihood under the parameter restriction
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- t0^kvek / factorial(kvek)
    for (j in 1:d){
      fvek <- fvek * t[j]^(x[i,j]-kvek) / factorial(x[i,j]-kvek)
    }
    fvek <- fvek * exp(-t0-sum(t))
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum)
}

# -log likelihood for a truncated Poisson sample, truncated at A.
dpolikut <- function(lam,x,A){
  n <- length(x)
  f <- -sum(log(lam ^ x / factorial(x) * exp(-lam))) + n * log(ppois(A, lam))
  return(f)
}

# -log likelihood for dependent Poisson distribution truncated at A
dpolikgqt <- function(t,A,x){
  n <- nrow(x)
  d <- ncol(x)
  liksum <- 0
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- t[1]^kvek / factorial(kvek)
    for (j in 1:d){
      fvek <- fvek * t[j]^(x[i,j]-kvek) / factorial(x[i,j] - kvek)
    }
    liksum <- liksum + log(sum(fvek))
  }
  tpsum <- 0
  for(u in 0:A){
    zvek <- c(0, cumsum(rep(1, A-u)))
    term <- t[1]^u / factorial(u)
    for(j in 1:d){
      term <- term * sum(t[j]^zvek / factorial(zvek))
    }
    tpsum <- tpsum + term
  }
  f <- -liksum + n * log(tpsum)
  return(f)
}

# -log likelihood for dependent Zero Inflated Poisson (ZIP) in dimension d
dziplikg <- function(plam,x){
  x <- as.matrix(x)
  plam <- as.vector(plam)
  n <- nrow(x)
  d <- ncol(x)
  d1 <- (d+1)*2
  liksum <- 0
  plammat <- matrix(plam, nrow = (d+1), ncol = 2, byrow = F)
  p <- plammat[,1]   #Extracts the pi:s."
  lam <- plammat[,2] #Extracts the lambdas."
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- zippdf(kvek, p[1], lam[1])
    #zippdf gives the pdf of zip for a vector of x:es"
    for (j in 1:d){
      fvek <- fvek*zippdf(x[i,j] - kvek, p[j+1], lam[j+1])
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum )
}

# -log likelihood for a dependent truncated Zero Inflated Poisson (ZIP) in dimension d
dziptlikg <- function(plam,x,A){
  x <- as.matrix(x)
  plam <- as.vector(plam)
  n <- nrow(x)
  d <- ncol(x)
  liksum <- 0
  plammat <- matrix(plam, nrow = (d+1), ncol = 2, byrow = F)
  p <- plammat[,1]     # Extracts the pi:s."
  lam <- plammat[,2]   # Extracts the lambdas."
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- ziptpdf(kvek, p[1], lam[1], A)
    for (j in 1:d){
      fvek <- fvek * ziptpdf(x[i,j] - kvek, p[j+1], lam[j+1], A)
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum )
}

# The pdf for Zero Inflated Poisson (ZIP) Distribution
zippdf <- function(x,p,lam){
  f <- p*(x==0) + (1-p)*exp(-lam)*lam^x/factorial(x)
  return(f)
}

# -log likelihood for univariate Zero Inflated Poisson (ZIP) Distribution
zippdf1 <- function(plam,x){
  p <- plam[1]
  lam <- plam[2]
  fvek <- p*(x == 0) + (1-p)*exp(-lam)*lam^x / factorial(x)
  return(-sum(log(fvek)))
}

# The pdf for Truncated Zero Inflated Poisson Distribution
ziptpdf <- function(x,p,lam,A){
  kvek <- c(0:A)
  pA <- p+(1-p)*exp(-lam)*sum(lam^kvek/factorial(kvek)) # truncation probability
  f <- (p*(x==0)+(1-p)*exp(-lam)*lam^x/factorial(x))/pA
  return(f)
}

# -log likelihood for univariate Truncated Zero Inflated Poisson Distribution
ziptpdf1 <- function(plam,x,A){
  p <- plam[1]
  lam <- plam[2]
  kvek <- c(0:A)
  pA <- p+(1-p)*exp(-lam)*sum(lam^kvek/factorial(kvek))    # truncation probability
  fvek <- (p*(x==0)+(1-p)*exp(-lam)*lam^x/factorial(x))/pA
  f <- sum(-log(fvek))
  return(f)
}

## Negative Binomial Distribution and Derivatives ####
# -log likelihood for dependent Negative Binomial in almant number of dimensions.
dnblikg1 <- function(rp,x){
  # Data x is n * d. Provides solutions with a wait value equal to the mean value for them
  # observed the variables. (These may be the right solutions, though
  #                          this must be displayed!) The parameter vector rp starts
  # with r0 (the r-parameter of the factor)
  # and below that come all the p-parameters. The dimension thus becomes d + 2;
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  d1 <- d + 2
  liksum <- 0
  r0 <- rp[1]
  p0 <- rp[2]
  p <- rp[3:d1]
  mx <- colMeans(x)
  r <- p/(1-p)*(mx-r0*(1-p0)/p0)
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- dnbinom(kvek, r0, p0)
    for (j in 1:d){
      fvek <- fvek*dnbinom(x[i,j] - kvek, r[j], p[j])
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum)
}

# -log likelihood for a sample from negative binomial.
dnbliku1 <- function(p,x){
  # Parametern r ges av m=r(1-p)/p.
  n <- length(x)
  m <- mean(x)
  r <- p/(1-p)*m
  liksum <- sum(log(dnbinom(x,r,p)))
  return(-liksum)
}

# -log likelihood for a truncated nagative binomial sample, truncated at A.
dnblikut <- function(rp,x,A){
  x <- as.matrix(x)
  n <- nrow(x)
  r <- rp[1]
  p <- rp[2]
  pA <- pnbinom(A,r,p)
  liksum <- 0
  fvek <- dnbinom(x,r,p)/pA
  liksum <- sum(log(fvek))
  return(-liksum)
}

# -log likelihood for dependent Negative Binomial in a general number of dimensions.
# Truncated at t.
dnblikgt <- function(rp,x,t){
  # The parameter vector rp is (d+1)*2, and data x is n*d. r is first in the vector, then p.
  n <- nrow(x)
  d <- ncol(x)
  d1 <- (d+1)*2
  rpmat <- matrix(rp, d+1, 2)
  r <- matrix(rpmat[, 1], ncol = 1)
  p <- matrix(rpmat[, 2], ncol = 1)
  liksum <- 0
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- (dnbinom(kvek, size = r[1, 1], prob = p[1, 1])) / pnbinom(t, size = r[1, 1], prob = p[1, 1])
    for (j in 1:d){
      fvek <- fvek * dnbinom(x[i,j] - kvek, r[j+1,1], p[j+1, 1]) / pnbinom(t, size = r[j+1, 1], prob = p[j+1, 1])
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum)
}

# The pdf for Zero Inflated Negative Binomial (ZINB).
zinbpdf <- function(x,pi,r,p){
  # Data is in a column vector x.
  # The parameters pi, r and p are scalars. (x==0) below is one for the entries
  # that are zero and zero otherwise.
  return(pi*(x==0) + (1-pi)*dnbinom(x,r,p))
}

# -log likelihood for one sample from zero inflated negative binomial.
dzinbliku <- function(pirp,x){
  n <- length(x)
  pi <- pirp[1]
  r <- pirp[2]
  p <- pirp[3]
  return(-sum(log(pi*(x == 0)+(1-pi)*dnbinom(x,r,p))))
}

# -log likelihood for dependent Zero Inflated Negative Binomial (ZINB) in dimension d.
dzinblikg <- function(pirp,x){
  # The parameter vector pirp is (d+1)*3, and data x is n*d. The pi:s are the
  # first entries in the vector (staring with the factor pi), then comes the
  # r in corresponding order, and then p.
  n <- nrow(x)
  d <- ncol(x)
  d1 <- (d+1)*3
  liksum <- 0
  pirpmat <- matrix(pirp, nrow = d+1, ncol = 3, byrow = F)
  pi <- pirpmat[,1]    # Extracts the pi:s.
  r <- pirpmat[,2]     # Extracts the r:s.
  p <- pirpmat[,3]     # Extracts the p:s.
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- zinbpdf(kvek, pi[1], r[1], p[1])
    for (j in 1:d){
      fvek <- fvek * zinbpdf(x[i,j] - kvek, pi[j+1], r[j+1], p[j+1])
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum)
}

# The pdf for Zero Inflated Negative Binomial (ZINB), truncated at A.
zinbtpdf <- function(x,pi,r,p,A){
  #Data is in a column vector x.
  #The parameters pi, r and p are scalars. (x==0) below is one for the entries
  #that are zero and zero otherwise.
  kvek <- c(0:A)
  pA <- sum(zinbpdf(kvek,pi,r,p))
  f <- (pi*(x==0)+(1-pi)*dnbinom(x,r,p))/pA
  return(f)
}

# -log likelihood for truncated zero inflated negative binomial, truncated at A.
dzinblikut <- function(pirp,x,A){
  n <- nrow(x)
  pi <- pirp[1]
  r <- pirp[2]
  p <- pirp[3]
  liksum <- 0
  kvek <- c(0:A)
  pA <- sum(zinbpdf(kvek,pi,r,p))
  liksum <- sum(log((pi*(x==0)+(1-pi)*dnbinom(x,r,p))/pA))
  return(-liksum)
}

# -log likelihood for dependent truncated Zero Inflated Negative Binomial (ZINB) in dimension d, truncated at A.
dzinblikgt <- function(pirp,x,A){
  # The parameter vector pirp is (d+1)*3, and data x is n*d. The pi:s are the
  # first entries in the vector (staring with the factor pi), then comes the
  # r in corresponding order, and then p.
  n <- nrow(x)
  d <- ncol(x)
  d1 <- (d+1)*3
  pirpmat <- matrix(pirp, nrow = d+1, ncol = 3)
  pi <- pirpmat[,1] # Extracts the pi's.
  r <- pirpmat[,2]  # Extracts the r's.
  p <- pirpmat[,3]  # Extracts the p's.
  liksum <- 0
  for (i in 1:n){
    m <- min(x[i,])
    kvek <- c(0, cumsum(rep(1,m)))
    fvek <- zinbtpdf(kvek, pi[1], r[1], p[1], A)
    for (j in 1:d){
      fvek <- fvek * zinbtpdf(x[i,j]-kvek, pi[j+1], r[j+1], p[j+1], A)
    }
    liksum <- liksum + log(sum(fvek))
  }
  return(-liksum)
}

