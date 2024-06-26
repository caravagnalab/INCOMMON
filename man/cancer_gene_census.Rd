% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cancer_gene_census}
\alias{cancer_gene_census}
\title{Gene roles from COSMIC Cancer Gene Census}
\format{
A data frame with 733 rows and 2 columns:
\describe{
\item{gene}{Name of the gene (Hugo Symbol)}
\item{gene_role}{Tumour Suppressor Gene (TSG) or oncogene}
}
}
\source{
COSMIC Cancer Gene Census: https://cancer.sanger.ac.uk/census
}
\usage{
data(cancer_gene_census)
}
\description{
This table contains gene roles in cancer, as obtained from the COSMIC Cancer
Gene Census v98. Data were curated such that each gene in list is assigned
either TSG or oncogene.
}
\examples{
data(cancer_gene_census)
cancer_gene_census
}
\keyword{datasets}
