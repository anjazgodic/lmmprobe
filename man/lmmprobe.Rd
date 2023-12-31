% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmmprobe.R
\name{lmmprobe}
\alias{lmmprobe}
\title{Sparse high-dimensional linear mixed modeling with PaRtitiOned empirical Bayes ECM (LMM-PROBE) algorithm.}
\arguments{
\item{Y}{A matrix containing the outcome \code{Y}.}

\item{Z}{A matrix containing the sparse fixed-effect predictors on which to apply the lmmprobe algorithm.}

\item{V}{A matrix containing non-sparse predictors for the random effects. This matrix is currently only programmed for two scenarios. Scenario 1: only a random intercept, where V is a matrix with one column containing ID's and each unit has the same number of observations. Scenario 2: a random intercept and a random slope, where V is a matrix with two columns. The first column is ID and the second column is a continuous variable (e.g. time) for which a random slope is to be estimated. Each unit has the same number of observations.}

\item{alpha}{Type I error; significance level.}

\item{ep}{Value against which to compare convergence criterion, we recommend 0.05.}

\item{B}{The number of groups to categorize estimated coefficients in to calculate correlation \eqn{\rho}. We recommend five.}

\item{adj}{A factor multiplying Silverman’s 'rule of thumb' in determining the bandwidth for density estimation, same as the 'adjust' argument of R's density function. Default is three.}

\item{maxit}{Maximum number of iterations the algorithm will run for. Default is 10000.}

\item{cpus}{The number of CPUS user would like to use for parallel computations. Default is four.}

\item{LR}{A learning rate parameter \code{r}. Using zero corresponds to the implementation described in Zgodic et al.}

\item{C}{A learning rate parameter \code{c}. Using one corresponds to the implementation described in Zgodic et al.}
}
\value{
A list of the output of the lmmprobe function, including 

\code{beta_hat, beta_hat_var} MAP estimates of the posterior expectation (beta_hat) and variance (beta_hat_var) of the prior mean (\eqn{\beta}) of the regression coefficients assuming \eqn{\gamma=1}, 

\code{gamma} the posterior expectation of the latent \eqn{\gamma} variables, 

\code{preds} predictions of \eqn{Y},

\code{PI_lower, PI_upper} lower and upper prediction intervals for the predictions,

\code{sigma2_est} MAP estimate of the residual variance, 

\code{random_var} MAP estimate of the random effect(s) variance, 

\code{random_intercept} estimated random intercept terms,  

\code{random_slope} estimated random slope terms, if applicable.
}
\description{
Sparse high-dimensional linear mixed modeling with PaRtitiOned empirical Bayes ECM (LMM-PROBE) algorithm. Currently, the package offers functionality for two scenarios. Scenario 1: only a random intercept, each unit has the same number of observations; Scenario 2: a random intercept and a random slope, each unit has the same number of observations. We are actively expanding the package for more flexibility and scenarios.
}
\examples{
library(lmmprobe)
data(SLE)
ep <- 0.05
alpha <- 0.05
Y = SLE$Y
Z = SLE$Z
V = SLE$V
full_res <- lmmprobe(Y = Y, Z = Z, V = V, ep = ep, alpha = alpha)
}
\references{
Zgodic, A. and McLain, A. C. (2023). Sparse high-dimensional linear mixed modeling with a partitioned empirical Bayes ECM algorithm. arXiv preprint arXiv:2310.12285.
}
