% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_jags_package.R
\docType{package}
\name{jagstargets-package}
\alias{jagstargets-package}
\alias{jagstargets}
\title{jagstargets: Targets for JAGS Workflows}
\description{
Bayesian data analysis usually incurs long runtimes
and cumbersome custom code. A pipeline toolkit tailored to
Bayesian statisticians, the \code{jagstargets} R package leverages
\code{targets} and \code{R2jags} to ease this burden.
\code{jagstargets} makes it super easy to set up scalable
JAGS pipelines that automatically parallelize the computation
and skip expensive steps when the results are already up to date.
Minimal custom code is required, and there is no need to manually
configure branching, so usage is much easier than \code{targets} alone.
}
\seealso{
\url{https://docs.ropensci.org/jagstargets/}, \code{\link[=tar_jags]{tar_jags()}}
}
