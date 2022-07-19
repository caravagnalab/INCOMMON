#' Prepare input for analyses with \code{'TAPACLOTH'}.
#' @description
#' Process input data into an object of class \code{'TAPACLOTH'}, ready for
#' analyses with `run_classifier` and `estimate_purity` functions.
#' @param mutations a tibble with columns chromosome `chr`, start position `from`, end position `to`,
#'   reference `ref` and alternative `alt` alleles, coverage `DP`, number
#'   of reads with variant `NV`, variant allelic frequency `VAF` gene name `gene` 
#'   as Hugo Symbol.
#' @param sample Sample name.
#' @param purity Sample purity.
#' @return An object of class \code{TAPACLOTH}.
#' @export
init = function(mutations, sample, purity, gene_roles = cancer_gene_census){
  mutations = left_join(mutations, gene_roles, by = "gene")
  out = list(
    data = mutations,
    sample = sample,
    purity = purity
  )
  class(out) = "TAPACLOTH"
  return(out)
}
