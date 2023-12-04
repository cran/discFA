#' @title Discrete factor analysis with the truncated zero inflated Poisson distribution
#'
#' @param y Data, an n by d numeric matrix
#' @param A truncation point (Note that if the data is in Likert scale
#' starting from 1, then you should subtract 1 from the data and then use the
#' proposed negative binomial models.
#'
#' @return A list with entries
#' \item{AIC}{AIC value for the optimal model}
#' \item{indexmat}{Factors and variables in each factor}
#' \item{estpilam}{Estimated zero-inflated parameters for  for each factor}
#' \item{estlam}{Estimated parameters for each factor}
#' \item{estpimu}{Estimated zero-inflated parameters for each variable within each factor}
#' \item{estmu}{Estimated parameters for each variable within each factor}
#'
#' @importFrom methods as
#' @importFrom stats nlminb
#' @importFrom stats ppois
#' @importFrom stats cor dnbinom pnbinom var
#' @export
#'
#' @examples
#' dfzipt(zinb100_Data[1:50,], A = 6)
dfzipt <- function(y,A=NULL){
  y <- as.matrix(y)
  # negative correlation
  Cor.Mat <- cor(y,method="kendall")
  if(any(Cor.Mat < -.5)) message("There are some highly negative correlations among some variables so the findings may not be stable.")
  ## warning and stop messages
  if(is.null(A)) stop("Please provide a value for A, truncation point.")
  if (any(is.na(y))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data and run dfziptgs again.")
  if (!is.numeric(y)) stop("Factor analysis applies only to numerical variables.")
  if (any(y<0)) stop("There are negative values in the data, please treat them and run factor analysis again.")
  if(is.null(colnames(y))){
    colnames(y) <- paste0("Var", c(1:ncol(y)), sep ="")
  }
  cn <- colnames(y)
  n <- nrow(y)
  d <- ncol(y)
  if (d < 3) stop("Factor analysis requires at least three variables.")

  cl <- match.call()
  start_time <- Sys.time()
  stop <- 0

  m <- colMeans(y)
  indexmat <- matrix(0, d, d)
  indexmat[1,] <- 1:d     # Initial grouping, the (1,1,...,1) model.

  optlik <- rep(0, d)
  for (j in 1:d){
    eval_f1 <-  function(tt){
      return(ff = ziptpdf1(tt,y[,j],A))
    }
    optlik[j] <- nlminb(c(0.0, m[j]), eval_f1,
                        lower = c(rep(0,2)),
                        upper = c(1, Inf),
                        control = list(trace = FALSE, x.tol = 1e-6))$objective
  }
  lik1 <- sum(optlik)
  AIC1 <- 2/n * (lik1 + 2*d)

  likmat <- 1000000*matrix(1, d, d)
  antilik <- matrix(0,d,d)

  for (i in 1:(d-1)){
    for (j in (i+1):d){
      ytemp <- cbind(y[,i], y[,j])
      mtemp <- c(m[i], m[j])
      eval_f2 <-  function(tt){
        return(ff = dziptlikg(tt, ytemp, A))
      }
      likmat[i,j] <- nlminb(c(rep(1e-06,4), mtemp), eval_f2,
                            lower = c(rep(0,6)),
                            upper = c(rep(1,3), rep(Inf,3)))$objective
      antilik[i,j] <- optlik[i] + optlik[j]
    }
  }
  lik2 <- likmat + lik1 - antilik
  min1 <- apply(lik2, 2, min)
  min_lik2 <- arrayInd(which.min(lik2), dim(lik2))
  I <- min_lik2[1] ; J <- min_lik2[2]
  minlik <- lik2[I,J]
  AIC2 <- 2/n*(minlik+2*(d+1))

  if (AIC1 < AIC2){ #This means that the independence model (1,1,...,1) is the best.
    AICmin <- AIC1
    stop <- 1
  } else {
    indexmat1 <- indexmat
    optlik1 <- optlik
    indexmat1[2,I] <- J
    optlik1[I] <- likmat[I,J]

    if(J < d){
      for(j in (J+1):d){  # Moving down...
        optlik1[j-1] <- optlik[j]
        indexmat1[1,(j-1)] <- j
      }
    }
    indexmat1[1, d] <- 0
    indexmat <- indexmat1
    posvar <- colSums(indexmat > 0)
    #This gives a vector with the number of variables in the d different columns of indexmat."
    optlik <- optlik1
    lik1 <- sum(optlik) # Optimal minus log likelihood.
    AIC1 <- AIC2        # New optimal AIC."
  }
  step <- 2
  while(stop == 0 && step < d)   {
    # Minimize minus loglik for all possible jonings of groups of variables from the previous step.
    d1 <- d - step + 1
    likmat <- 1000000*matrix(1, d1, d1)
    antilik <- matrix(0, d1, d1)

    for (i in 1:(d1-1)){
      for (j in (i+1):d1){
        ytemp <- cbind(y[,indexmat[c(1:posvar[i]),i]], y[,indexmat[c(1:posvar[j]),j]])
        mtemp <- colMeans(ytemp)
        nvar <- posvar[i] + posvar[j] + 1
        eval_f3 = function(tt){
          return(ff = dziptlikg(tt,ytemp,A))
        }
        pars0 <- c(rep(1e-06, nvar+1), mtemp)
        likmat[i,j] <- nlminb(pars0, eval_f3,
                              lower = rep(0, length(pars0)),
                              upper = c(rep(1, nvar), rep(Inf, nvar)))$objective

        antilik[i,j] <- optlik[i] + optlik[j] - 2 *(1 - (posvar[i] > 1) - (posvar[j] > 1))
      } #inner for loop ends...
    } #outer for loop ends...

    lik2 <- likmat + minlik - antilik
    # minlik is the optimal minus log likelihood (corrected for the number of parameters) in the previous step.

    min1 <- apply(lik2, 2, min)
    min_lik2 <- arrayInd(which.min(lik2), dim(lik2))
    I <- min_lik2[1] ; J <- min_lik2[2]
    minlik <- lik2[I,J]
    ################################################
    posI <- posvar[I]
    posJ <- posvar[J]
    indexmat1 <- indexmat
    optlik1 <- optlik
    indexmat1[(posI+1):(posI+posJ), I] <- indexmat[(1:posJ), J]
    optlik1[I] <- likmat[I,J]
    if(J < d1){
      for(j in (J+1):d1){  # Moving down...
        optlik1[j-1] <- optlik[j]
        indexmat1[,(j-1)] <- indexmat[,j]
      }
    }
    indexmat1[,d1] <- 0
    posvar1 <- colSums(indexmat1 > 0)
    AIC2 <- (2/n)*(minlik + 2*(d + 1))
    if (AIC1 < AIC2){ # The previous model was better...
      AICmin <- AIC1
      stop <- 1
    } else{
      optlik <- optlik1
      posvar <- posvar1
      indexmat <- indexmat1
      AIC1 <- AIC2 # New optimal AIC.
      step <- step+1
    } # else ends..
  } # while ends...
  AICmin <- AIC1 # Gives the output AIC.
  j <- 1
  estpilam <- estlam <- c(rep(0, d))
  estpimu <- estmu <- matrix(0, nrow = d, ncol = d)
  stop <- 0
  while(stop == 0) {                     # indicates a non empty submodel
    nvar <- sum(indexmat[, j] > 0)       # number of variables in submodel
    vartemp <- indexmat[1:nvar, j]       # variable indeces in submodel
    ytemp <- as.matrix(y[, vartemp])     # corresponding observations
    mtemp <- colMeans(ytemp)
    nvar1 <- nvar + 1

    if(nvar == 1){
      ftemp1 <-  function(tt1){
        return(ff = ziptpdf1(tt1, ytemp, A))
      }

      OPT <- nlminb(c(1e-06, m[j]), ftemp1,
                    lower = c(0,0),
                    upper = c(1,Inf))
      est <- OPT$par
      estpimu[1,j] <- est[1]
      estmu[1,j] <- est[2]
    } else {
      ftemp2 <-  function(tt){
        return(ff = dziptlikg(tt, ytemp, A))
      }
      OPT2 <- nlminb(c(rep(1e-06, nvar+2), mtemp), ftemp2,
                     lower = c(rep(0, 2*nvar+2)),
                     upper = c(rep(1, nvar1),rep(Inf, nvar1)))

      est <- OPT2$par
      estpilam[j] <- est[1]
      estlam[j] <- est[nvar+2]
      estpimu[1:nvar, j] <- est[2:nvar1]
      estmu[1:nvar, j] <- est[(nvar+3):(2*nvar1)] # estimated mu's for jth submodel
    } # end of else

    j <- j + 1
    if(j > d){
      stop <- 1
    } else if (indexmat[1,j] == 0) {
      stop <- 1
    }
  } # end of while loop
  finish_time <- Sys.time()
  total_time <- finish_time-start_time

  not_zero <- as.vector(colSums(as.matrix(indexmat!=0)))
  not_zero_index <- which(not_zero!=0)
  not_zero <- not_zero[c(not_zero_index)]

  if(sum(colSums(indexmat[-1,]))==0) {
    row.zeros <- 1
    n.all.zeros <- d
  } else if (which(colSums(indexmat) == 0)[1]==2){
    row.zeros <- d
    n.all.zeros <- which(colSums(indexmat==0) == d)[1]-1
  } else {
    row.zeros <- which(rowSums(indexmat) == 0)[1]-1
    n.all.zeros <- which(colSums(indexmat==0) == d)[1]-1
  }

  indexmat.final <- matrix(indexmat[1:row.zeros,1:n.all.zeros], ncol=n.all.zeros)

  indexmat.final[which(indexmat.final!=0)] <- cn[as.vector(indexmat.final)]
  colnames(indexmat.final) <- paste0("Factor", c(1:ncol(indexmat.final)), sep ="")
  rownames(indexmat.final) <- paste0(c(1:nrow(indexmat.final)), sep ="")

  estmu.final <- matrix(estmu[1:row.zeros,1:n.all.zeros], ncol=n.all.zeros)
  colnames(estmu.final) <- paste0("Factor", c(1:ncol(estmu.final)), sep ="")
  rownames(estmu.final) <- paste0(c(1:nrow(estmu.final)), sep ="")

  estpimu.final <- matrix(estpimu[1:row.zeros,1:n.all.zeros], ncol=n.all.zeros)
  colnames(estpimu.final) <- paste0("Factor", c(1:ncol(estpimu.final)), sep ="")
  rownames(estpimu.final) <- paste0(c(1:nrow(estpimu.final)), sep ="")

  result <- list(AIC=AICmin,
                 indexmat=indexmat.final,
                 estpimu=estpimu.final,
                 estmu=estmu.final,
                 estpilam=estpilam,
                 estlam=estlam,timing=total_time,n=n, d=d)
  result$call <- cl
  result$model <- not_zero
  class(result) <- "dfzipt"
  result
} # dfzipt ends...
#' @export
print.dfzipt <- function(x, digits = 4, ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

  if(setequal(x$model, rep(1,x$d))){message("Independent model!\n")}
  cat("This is a (",toString(paste0(x$model)),") model.\n",sep="")

  if(x$AIC == Inf){
    message("The AIC is Inf and discrete factor analysis may not be feasible for this data.\n")
  } else {
    cat("\nAIC value is ", x$AIC, ".\n", sep="")
  }

  cat("\nFactors and variables in each factor:\n")
  print(ifelse(x$indexmat == 0, "", x$indexmat), quote = FALSE, ...)

  cat("\nEstimated zero-inflated parameters for each variable within each factor:\n")
  print(ifelse(x$indexmat == 0, "", round(x$estpimu, digits)), quote = FALSE, ...)
  #print(round(x$estpimu, digits), quote = FALSE, ...)

  cat("\nEstimated parameters for each variable within each factor:\n")
  print(ifelse(x$indexmat == 0, "", round(x$estmu, digits)), quote = FALSE, ...)
  #print(round(x$estmu, digits), quote = FALSE, ...)

  cat("\nEstimated zero-inflated parameters for each factor:\n")
  cat(round(x$estpilam[which(x$estpilam!=0)], digits=digits))
  #cat(round(x$estpilam, digits=digits))

  cat("\n\nEstimated parameters for factors:\n")
  cat(round(x$estlam[which(x$estlam!=0)], digits=digits))
  #cat(round(x$estlam, digits=digits))

  cat("\n\nTiming:\n")
  print(x$timing, digits = digits, quote = FALSE, ...)
  invisible(x)
}
