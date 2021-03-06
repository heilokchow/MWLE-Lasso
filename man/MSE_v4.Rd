% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Collection of Functions.R
\name{MSE_v4}
\alias{MSE_v4}
\title{MSE of changed ability follows step function with two teams involved and no home effect}
\usage{
MSE_v4(a, n, a1, k1, p, m)
}
\arguments{
\item{a}{Exponetial decay rate.}

\item{n}{Number of matched.}

\item{a1}{The ability score at the beginning.}

\item{k1}{Gap of step changed ability}

\item{p}{Lasso penalty level = s/max(s)}

\item{m}{Proportion of most recent unchanged ability}
}
\value{
Bias, Variance and MSE
}
\description{
MSE of changed ability follows step function with two teams involved and no home effect
}
