#' @title Discrete factor analysis with the Poisson distribution
#'
#' @param y Data, an n by d numeric matrix
#'
#' @return A list with entries
#' \item{AIC}{AIC value for the optimal model}
#' \item{indexmat}{Factors and variables in each factor}
#' \item{estlambda}{Estimated parameters for factors}
#' \item{estmu}{Estimated parameters for each variable within each factor}
#'
#' @export
#' @importFrom methods as
#' @importFrom stats nlminb
#' @importFrom stats ppois
#' @importFrom stats cor dnbinom pnbinom var
#'
#' @examples
#' dfp(car_data[,2:9])
dfp <- function(y){
  y <- as.matrix(y)
  # negative correlation
  Cor.Mat <- cor(y,method="kendall")
  if(any(Cor.Mat < -.5)) message("There are some highly negative correlations among some variables so the findings may not be stable.")
  ## warning and stop messages
  if (any(is.na(y))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data and run dfpogs again.")
  if (!is.numeric(y)) stop("Factor analysis applies only to numerical variables.")
  if (any(y<0)) stop("There are negative values in the data, please treat them and run factor analysis again.")
  if(is.null(colnames(y))){colnames(y) <- paste0("Var", c(1:ncol(y)), sep ="")}
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
    optlik[j] <- dpoliku(y[,j]) # Gives univariate maximum likelihood.
  }
  lik1 <- sum(optlik)
  AIC1 <- (2/n) * (lik1 + d)

  likmat <- 1000000*matrix(1, d, d) # To prevent that no element with i>j is choosen in the minimization.
  antilik <- matrix(0,d,d)
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      ytemp <- cbind(y[,i], y[,j])
      eval_f1 = function(t0){
        return(ff = dpolikgq(t0, ytemp))
      }
      likmat[i,j] <- nlminb(0, eval_f1,
                            lower = 0, upper = Inf,
                            control = list(trace = FALSE, rel.tol = 1e-6))$objective
      antilik[i,j] <- optlik[i] + optlik[j]
    }
  }
  lik2 <- likmat + lik1 - antilik

  # Below: All optimal minus logliks are saved in a d vector.
  # Only join the variables/groups that need to be joined.
  # Move the optimal minus logliks and indices in indexmat down.
  # This way, we don not need to redo any maximizations already done in the
  # next step.
  min1 <- apply(lik2, 2, min)
  min_lik2 <- arrayInd(which.min(lik2), dim(lik2))
  I <- min_lik2[1] ; J <- min_lik2[2]
  minlik <- lik2[I,J]
  AIC2 <- 2/n*(minlik+d+1)

  if (AIC1 < AIC2){ # This means that the independence model (1,1,...,1) is the best.
    AICmin <- AIC1
    stop <- 1
  } else{
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
    indexmat1[1,d] <- 0
    indexmat <- indexmat1
    posvar <- colSums(indexmat>0)

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
        eval_f2 <- function(t0){
          return(ff = dpolikgq(t0,(ytemp)))
        }
        likmat[i,j] <- nlminb(0, eval_f2,
                              lower = 0, upper = Inf,
                              control = list(trace = FALSE, abs.tol = 1e-6))$objective

        antilik[i,j] <- optlik[i] + optlik[j] - 1 + (posvar[i] > 1) + (posvar[j] > 1)
      } #inner for loop ends...
    } #outer for loop ends...

    lik2 <- likmat + minlik - antilik
    # minlik is the optimal minus log likelihood (corrected for the number of parameters) in the previous step.

    min1 <- apply(lik2, 2, min)
    min_lik2 <- arrayInd(which.min(lik2), dim(lik2))
    I <- min_lik2[1] ; J <- min_lik2[2]
    minlik <- lik2[I,J]

    posI <- posvar[I]
    posJ <- posvar[J]
    indexmat1 <- indexmat
    optlik1 <- optlik
    indexmat1[(posI+1):(posI+posJ),I] <- indexmat[(1:posJ),J]
    optlik1[I] <- likmat[I,J]

    if(J < d1){
      for(j in (J+1):d1){  # Moving down...
        optlik1[j-1] <- optlik[j]
        indexmat1[,(j-1)] <- indexmat[,j]
      }
    }

    indexmat1[,d1] <- 0
    posvar1 <- colSums(indexmat1 > 0)
    AIC2 <- 2/n*(minlik+d+1) # OBS: minlik is already corrected for the changed number of parameters compared to the previous step.

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
  estlambda <- c(rep(0,d))
  estmu <- matrix(0, nrow = d, ncol = d)
  stop <- 0
  while(stop == 0) {                     # indicates a non empty submodel
    nvar <- sum(indexmat[, j] > 0)       # number of variables in submodel
    vartemp <- indexmat[1:nvar, j]       # variable indeces in submodel
    ytemp <- as.matrix(y[, vartemp])     # corresponding observations

    if(nvar == 1){
      liktemp <- dpoliku(ytemp)
      estmu[1,j] <- mean(ytemp)
    } else {
      eval_f3 = function(tt){
        return(ff = dpolikgq(tt, ytemp))
      }
      OPT <- nlminb(0, eval_f3,
                    lower = 0, upper = Inf,
                    control = list(trace = FALSE, abs.tol = 1e-6))

      estlambda[j] <- OPT$par
      estmu[1:nvar, j] <- colMeans(ytemp) - OPT$par
    }
    j <- j + 1
    if(j > d){
      stop <- 1
    } else if (indexmat[1,j] == 0) {
      stop <- 1
    }
  }
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

  result <- list(AIC=AICmin,indexmat=indexmat.final,
                 estlambda=estlambda,estmu=estmu.final,
                 timing=total_time, n=n, d=d)
  result$call <- cl
  result$model <- not_zero
  class(result) <- "dfp"
  result
} # dfp function ends...

#' @export
print.dfp <- function(x, digits = 4, ...){
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

  cat("\nEstimated parameters for each variable within each factor:\n")
  print(ifelse(x$estmu == 0, "", round(x$estmu, digits)), quote = FALSE, ...)

  cat("\nEstimated parameters for factors:\n")
  cat(round(x$estlambda[which(x$estlambda!=0)], digits=digits))
  #cat(round(x$estlambda, digits=digits))

  cat("\n\nTiming:\n")
  print(x$timing, digits = digits, quote = FALSE, ...)
  invisible(x)
}
