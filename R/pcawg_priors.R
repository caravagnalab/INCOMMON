#' @title PCAWG empirical priors.
#'
#' @description Prior distribution of copy number and mutation multiplicity from PCAWG.
#' 
#' @format A data frame with 4466 rows and 4 columns:
#' \describe{
#'   \item{gene}{Name of the gene (Hugo Symbol).}
#'   \item{tumor_type}{Tumor type.}
#'   \item{label}{INCOMMON class (<total number of copies>N (Mutated: <mutation multiplicity>N)).}
#'   \item{p}{Gene and tumor type specific prior probability.}
#' }
#' @source Validated copy number calls from PCAWG: https://doi.org/10.5281/zenodo.6410935
#' @examples
#' \dontrun{
#'   head(pcawg_priors)
#' }