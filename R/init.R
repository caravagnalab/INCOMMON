#' Prepare input for analyses with \code{'INCOMMON'}.
#' @description
#' Process input data into an object of class \code{'INCOMMON'}, ready for
#' dwonstream analyses (e.g. `classify_sample`).
#' @param mutations a tibble with columns mutant chromosome `chr`, mutation start position `from`, mutation end position `to`,
#'   reference allel `ref`, alternative allele `alt`, integer sequencing depth `DP`, integer number
#'   of reads with variant `NV`, variant allelic frequency `VAF` and gene name `gene`
#'   as Hugo Symbol.
#' @param sample Sample name or ID.
#' @param tumor_type Tumor type of the sample.
#' @param purity Purity of the sample.
#' @param gene_roles A tibble reporting a list of `gene` names and associated
#'  `gene_role` among "oncogene", "TSG" and "fusion". The default is taken from COSMIC Cancer Gene Census.
#' @return An object of class \code{INCOMMON}.
#' @export
#' @examples 
#' input = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity, tumor_type = example_data$tumor_type)
#' print(input)
init = function(mutations, sample, tumor_type, purity, gene_roles = INCOMMON::cancer_gene_census){
  mutations = dplyr::left_join(mutations, gene_roles, by = "gene")
  out = list(
    data = mutations,
    sample = sample,
    purity = purity,
    tumor_type = tumor_type
  )
  class(out) = "INCOMMON"
  return(out)
}
