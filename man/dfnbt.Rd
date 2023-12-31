% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfnbt.R
\name{dfnbt}
\alias{dfnbt}
\title{Discrete factor analysis for the truncated negative binomial distribution (with right truncation at A)}
\usage{
dfnbt(y, A)
}
\arguments{
\item{y}{Data, an n by d numeric matrix}

\item{A}{truncation point (Note that if the data is in Likert scale
starting from 1, then you should subtract 1 from the data and then use the
proposed negative binomial models.}
}
\value{
A list with entries
\item{AIC}{AIC value for the optimal model}
\item{indexmat}{Factors and variables in each factor}
\item{estr0}{Estimated value of r for the factor}
\item{estp0}{Estimated value of p for the factor}
\item{estr}{Estimated value of r for the observation}
\item{estp}{Estimated value of r for the observation}
}
\description{
Discrete factor analysis for the truncated negative binomial distribution (with right truncation at A)
}
\examples{
dfnbt(zinb100_Data[1:40,1:3], A = 6)
}
