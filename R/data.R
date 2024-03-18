#' @docType data
#'
#' @usage data(example_data)
#'
#' @description Genomic data from sample P-0002081-T01-IM3 of the MSK-MET cohort.
#' 
#' @format A named list with the following components:
#' \describe{
#'   \item{mutations}{A data frame containing mutations.}
#'   \item{sample}{A string indicating the name of the sample.}
#'   \item{tumor_type}{A string indicating the tumor type.}
#'   \item{purity}{A numeric indicating the fraction of tumor in sample.}
#' }
#'
#' The columns in the data frame (\code{mutations}) are:
#' \describe{
#'   \item{chr}{Chromosome harbouring the mutation.}
#'   \item{from}{Genomic start position of the mutation.}
#'   \item{to}{Genomic end position of the mutation.}
#'   \item{ref}{Reference allele.}
#'   \item{alt}{Alternative allele.}
#'   \item{gene}{Gene harbouring the mutation.}
#'   \item{NV}{Number of reads with the variant}
#'   \item{DP}{Total number of reads (sequencing depth)}
#'   \item{VAF}{Variant allele frequency of the mutation.}
#' }
#'
#' @details This list consists the following components:
#' \itemize{
#'   \item \code{mutations}: A data frame containing mutations from sample P-0002081-T01-IM3 
#'   of the MSK-MET cohort with columns `chr`, `from`, `to`, `ref`, `alt`, `gene`, `NV`, `DP` and `VAF`.
#'   \item \code{sample}: The name of the sample.
#'   \item \code{tumor_type}: The tumor type of the sample.
#'   \item \code{purity}: The purity of the tumor sample (fraction).
#' }
#' 
#' @examples
#' data(example_data)
#' example_data
#' \dontrun{
#' # Create an example data
#' example <- list(
#' mutations <- data.frame(
#'  chr = c("chr1", "chr14", "chr17"),  # Chromosome column
#'  from = c(16265908, 105246551, 7578503),  # Start position column
#'  to = c(16265908, 105246551, 7578518),  # End position column
#'  ref = c("A", "C", "CAGGGCAGGTCTTGGC"),  # Reference allele column
#'  alt = c("T", "T", "-"),  # Alternate allele column
#'  gene = c("SPEN", "AKT1", "TP53")  # Gene column
#' ),
#'   sample = "P-0000004-T01-IM3",
#'   purity = 0.5
#' )
#' }
#'
#' @source MSK-MET at cBioPortal https://www.cbioportal.org/study/summary?id=msk_met_2021

#' @docType data
#'
#' @usage data(pcawg_priors)
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
#' data(pcawg_priors)
#' pcawg_priors
 
#' @docType Data
#'
#' @description This table contains gene roles in cancer, as obtained from the COSMIC Cancer
#' Gene Census v98. Data were curated such that each gene in list is assigned 
#' either TSG or oncogene.
#' 
#' @usage data(cancer_gene_census)
#' 
#' @format A data frame with 733 rows and 2 columns:
#' \describe{
#'   \item{gene}{Name of the gene (Hugo Symbol)}
#'   \item{gene_role}{Tumour Suppressor Gene (TSG) or oncogene}
#' }
#' @source COSMIC Cancer Gene Census: https://cancer.sanger.ac.uk/census
#' @examples
#' data(cancer_gene_census)
#' cancer_gene_census