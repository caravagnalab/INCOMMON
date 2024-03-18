#' @title Gene roles (TSG/oncogene) from COSMIC Cancer Gene Census
#'
#' @description This table contains gene roles in cancer, as obtained from the COSMIC Cancer
#' Gene Census v98. Data were curated such that each gene in list is assigned 
#' either TSG or oncogene.
#' 
#' @format A data frame with 733 rows and 2 columns:
#' \describe{
#'   \item{gene}{Name of the gene (Hugo Symbol)}
#'   \item{gene_role}{Tumour Suppressor Gene (TSG) or oncogene}
#' }
#' @source COSMIC Cancer Gene Census: https://cancer.sanger.ac.uk/census
#' @examples
#' \dontrun{
#'   head(cancer_gene_census)
#' }