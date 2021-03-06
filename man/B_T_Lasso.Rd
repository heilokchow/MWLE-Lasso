% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Collection of Functions.R
\name{B_T_Lasso}
\alias{B_T_Lasso}
\title{Bradley-Terry Likelihood Function with Lasso Peanlty's optimization}
\usage{
B_T_Lasso(dataframe, lambda, wij)
}
\arguments{
\item{dataframe}{Dataframe with 4 columns. First column is the home teams
Second column is the away teams.
Third column is the number of wins of home teams.
Fourth column is the number of wins of away teams.}

\item{lambda}{Lasso penalty of no statistical mean, a larger choice of lambda means higher penalty.
Usually the best penalty is chosen from Cross-Validation, where in-model's prediction power is maximized.}

\item{wij}{The weights added to the Lasso Penalty. Can be manually setted or determined using
function B_T_Wij}
}
\description{
Bradley-Terry Likelihood Function with Lasso Peanlty's optimization
}
