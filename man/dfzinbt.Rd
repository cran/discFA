% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfzinbt.R
\name{dfzinbt}
\alias{dfzinbt}
\title{Discrete factor analysis with the zero inflated truncated negative binomial distribution.}
\usage{
dfzinbt(y, A, tol = 1e-06)
}
\arguments{
\item{y}{Data, an n by d numeric matrix}

\item{A}{truncation point (Note that if the data is in Likert scale
starting from 1, then you should subtract 1 from the data and then use the
proposed negative binomial models.}

\item{tol}{tolerance value for optimizations}
}
\value{
A list with entries.
\item{AIC}{AIC value for the optimal model}
\item{indexmat}{Factors and variables in each factor}
\item{estpi0}{Estimated value of p for the zero inflated part in the negative binomial distributed factor}
\item{estr0}{Estimated value of r the negative binomial distributed factor(s)}
\item{estp0}{Estimated value of p the negative binomial distributed factor(s)}
\item{estpi}{Estimated parameters for the zero inflated part in the negative binomial distributed observations(s)}
\item{estr}{Estimated value of r negative binomial distributed observations(s)}
\item{estp}{Estimated value of p negative binomial distributed observations(s)}
}
\description{
Discrete factor analysis with the zero inflated truncated negative binomial distribution.
}
\examples{
dfzinbt(zinb100_Data[1:20,1:3], A = 6)
}
