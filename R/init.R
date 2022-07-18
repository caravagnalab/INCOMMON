#' #' Prepare input for analyses with \code{'TAPACLOTH'}.
#' @description
#' Cast input data in the format requested by \code{'TAPACLOTH} `run_classifier`
#' and `estimate_purity` functions.
#' @param mutations A tibble listing mutations called in the sample
#' in the form of chromosome `chr`, start position `from`, end position `to`,
#' reference `ref` and alternative `alt` alleles, gene name `gene` as Hugo Symbol.
#' @param sample Sample name.
#' @param purity Sample purity.
#' @return A list with mutations including gene role as annotated in COSMIC
#' cancer gene census, sample name and purity, in the format requested for 
#' \code{TAPACLOTH} analyses.
#' @export
init = function(mutations, sample, purity){
  mutations = left_join(mutations, cancer_gene_census, by = "gene")
  out = list(
    data = mutations,
    sample = sample,
    purity = purity
  )
}