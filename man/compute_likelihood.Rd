% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_likelihood.R
\name{compute_likelihood}
\alias{compute_likelihood}
\title{Compute model likelihood}
\usage{
compute_likelihood(NV, DP, prob, rho)
}
\arguments{
\item{NV}{Number of reads with the variant.}

\item{DP}{Sequencing coverage of the mutated genome site.}

\item{prob}{Success probability (expected VAF).}

\item{rho}{The over-dispersion parameter.}
}
\value{
A vector of probability densities (from NV = 1 to NV = DP).
}
\description{
Compute model likelihood
}
\examples{
compute_likelihood(
NV = 170,
DP = 200,
prob = 0.5,
rho = 0.01
)
}
