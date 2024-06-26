% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox_fit.R
\name{cox_fit}
\alias{cox_fit}
\title{Fit multivariate Cox regression model based on INCOMMON classes.}
\usage{
cox_fit(
  x,
  gene,
  tumor_type,
  survival_time,
  survival_status,
  covariates = c("age", "sex", "tmb"),
  tmb_method = "median"
)
}
\arguments{
\item{x}{An object of class \code{'INCOMMON'} containing the classification results
as produced by function \code{classify}.}

\item{gene}{The gene on which patient's stratification is based.}

\item{tumor_type}{The tumor type of patients to stratify.}

\item{survival_time}{The variable in \code{clincal_data} to be used as survival time.}

\item{survival_status}{The variable in \code{clincal_data} to be used as survival status.}

\item{covariates}{The other covariates to be used in the mutlivariate regression.}

\item{tmb_method}{The method to define the reference value for tumor mutational burden TMB}
}
\value{
An object of class \code{'INCOMMON'} containing an additional object \code{survival}.
}
\description{
Fit multivariate Cox regression model based on INCOMMON classes.
}
\examples{
# First load example classified data
data(MSK_classified)
# Perform Cox regression based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
MSK_classified = cox_fit(x = MSK_classified, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS', covariates = c('age', 'sex', 'tmb'), tmb_method = ">10")
}
