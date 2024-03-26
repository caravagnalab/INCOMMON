#' Prepare input for analyses with \code{'INCOMMON'}.
#' @description
#' Process input data into an object of class \code{'INCOMMON'}, ready for
#' downstream analyses (e.g. `classify`).
#' @param genomic_data a data table of annotated mutations with columns sample name `sample`, mutant chromosome `chr`, mutation start position `from`, mutation end position `to`,
#'   reference allele `ref`, alternative allele `alt`, integer sequencing depth `DP`, integer number
#'   of reads with variant `NV`, variant allele frequency `VAF` and gene name `gene`
#'   as Hugo Symbol,  protein sequence of the variant in HGVS recommended format, preferably 1-letter amino-acid code `HGVSp_Short`.
#' @param clinical_data a data table of clinical data with compulsory matching sample names `sample` and sample purity `purity`, and optional
#' clinical features like tumor type ONCOTREE code `tumor_type` (required for tumor specific priors), overall survival status `OS_STATUS` and time `OS_MONTHS`
#' (required for survival analysis), `SAMPLE_TYPE` (Primary or Metastasis) and  number of metastases `MET_COUNT` (required for metastatic propensity analysis),
#' metastatic site `METASTATIC_SITE` (required for metastatic tropism analysis), plus any other useful covariate.
#' @param gene_roles A data table reporting `gene` names and associated `gene_role` ("oncogene" or "TSG"). The default is taken from COSMIC Cancer Gene Census v98.
#' @return An object of class \code{INCOMMON}.
#' @export
#' @examples
#' input = init(genomic_data = MSK_genomic_data, clinical_data = MSK_clinical_data)
#' print(input)
init = function(genomic_data,
                clinical_data,
                gene_roles = INCOMMON::cancer_gene_census) {
  # Add gene roles
  genomic_data = genomic_data %>%
    dplyr::left_join(gene_roles, by = "gene")
  # Create object as list
  out = list(genomic_data = genomic_data,
             clinical_data = clinical_data)
  # Add class specification
  class(out) = "INCOMMON"

  ## Print input data
  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::black("[ INCOMMON ] ")),
      'Inference of copy number and mutation multiplicity in oncology'
    )
  )

  # Genomic data summary

  n_samples = out$genomic_data$sample %>% unique() %>% length()
  n_mutations = out$genomic_data %>% nrow()
  n_genes = out$genomic_data$gene %>% unique() %>% length()

  na_gene_names = dplyr::filter(out$genomic_data, is.na(gene)) %>% nrow()
  na_role_genes = dplyr::filter(out$genomic_data, is.na(gene_role)) %>% dplyr::pull(gene) %>% unique() %>% length()
  no_reads_mutations = dplyr::filter(out$genomic_data, is.na(NV) | is.na(DP)) %>% nrow()
  no_reads_samples = dplyr::filter(out$genomic_data, is.na(NV) | is.na(DP)) %>% dplyr::pull(sample) %>% unique() %>% length()

  cli::cli_h2("Genomic data")

  cli::cli_alert_success(
    "Found {.field {n_samples}} samples, with {.field {n_mutations}} mutations in {.field {n_genes}} genes"
  )

  if(no_reads_mutations != 0)
    cli::cli_alert_warning("No read counts found for {.field {no_reads_mutations}} mutations in {.field {no_reads_samples}} samples")

  if(na_gene_names != 0)
    cli::cli_alert_warning("Gene name not provided for {.field {na_gene_names}} mutations")

  if(na_role_genes != 0)
    cli::cli_alert_warning("{.field {na_role_genes}} genes could not be assigned a role (TSG or oncogene)")

  mean_depth = genomic_data$DP %>% mean(na.rm = T) %>% as.integer()

  cli::cli_alert_info('Average sequencing depth: {.field {mean_depth}}')

  # Clinical data summary

  clinical_features = colnames(out$clinical_data)

  cli::cli_h2("Clinical data")

  cli::cli_alert_info(
    "Provided clinical features:"
  )
  cat('\n')

  missing_features = c()

  for(x in c('sample', 'purity')){
    if(x %in% clinical_features) {
      cli::cli_bullets(c('v' = paste(x, '(required for classification)')))
    } else{
      cli::cli_bullets(c('x' = paste(x, '(required for classification)')))
      missing_features = append(x, missing_features)
    }
  }

  for(x in setdiff(clinical_features, c('sample', 'purity'))){cli::cli_bullets(c('v' = paste(x)))}

  cat('\n')


  if(!is.null(missing_features)) cli::cli_alert_warning('Missing features {.field {missing_features}} required for classifcation')

  matching_samples = intersect(genomic_data$sample, clinical_data$sample) %>% unique() %>% length()
  unmatched_samples = c(setdiff(genomic_data$sample, clinical_data$sample), setdiff(clinical_data$sample, genomic_data$sample)) %>% unique() %>% length()

  cli::cli_alert_success("Found {.field {matching_samples}} matching samples")
  if(unmatched_samples != 0) cli::cli_alert_danger("Found {.field {unmatched_samples}} matching samples") else  cli::cli_alert_success("No mismatched samples")

  mean_purity = clinical_data$purity %>% mean(na.rm = T) %>% round(2)

  cli::cli_alert_info('Average sample purity: {.field {mean_purity}}')

  return(out)
}
