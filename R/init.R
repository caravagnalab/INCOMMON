#' Prepare input for analyses with \code{'INCOMMON'}.
#' @description
#' Process input data into an object of class \code{'INCOMMON'}, ready for
#' dwonstream analyses (e.g. `classify_sample`).
#' @param mutations a tibble with columns mutant chromosome `chr`, mutation start position `from`, mutation end position `to`,
#'   reference allel `ref`, alternative allele `alt`, integer sequencing depth `DP`, integer number
#'   of reads with variant `NV`, variant allelic frequency `VAF` and gene name `gene`
#'   as Hugo Symbol.
#' @param sample Sample name.
#' @param purity Sample purity.
#' @param gene_roles A tibble reporting a list of `gene` names and associated
#'  `gene_role` among "oncogene", "TSG" and "fusion". The default is taken from COSMIC Cancer Gene Census.
#' @return An object of class \code{TAPACLOTH}.
#' @export
#' @examples 
#' input = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity)
#' print(input)
init = function(mutations, sample, purity, gene_roles = TAPACLOTH::cancer_gene_census){
  mutations = dplyr::left_join(mutations, gene_roles, by = "gene")
  out = list(
    data = mutations,
    sample = sample,
    purity = purity
  )
  class(out) = "TAPACLOTH"
  return(out)
}
