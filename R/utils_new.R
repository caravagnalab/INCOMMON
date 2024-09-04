check_input = function(x){
  if(genomic_data(x)$sample %>% class() != "character")
    cli::cli_abort("Sample names must be characters, classification aborted.")

  if(genomic_data(x)$gene %>% class() != "character")
    cli::cli_abort("Gene names must be characters, classification aborted.")

  if(genomic_data(x)$gene_role %>% class() != "character" |
     setdiff(genomic_data(x)$gene_role, c("TSG", "oncogene", NA)) %>% length() != 0)
    cli::cli_abort("Gene roles must be \'TSG\' or \'onocgene\', classification aborted.")

  if(genomic_data(x)$chr %>% class() != "character")
    cli::cli_abort("Chromosome names must be characters, classification aborted.")

  if(genomic_data(x)$ref %>% class() != "character")
    cli::cli_abort("Reference/Alternative allele names must be characters, classification aborted.")

  if(genomic_data(x)$alt %>% class() != "character")
    cli::cli_abort("Reference/Alternative allele names must be characters, classification aborted.")

  if(genomic_data(x)$VAF %>% class() != "numeric")
    cli::cli_abort("VAF must be numeric, classification aborted.")

  if(genomic_data(x)$DP %>% class() != "integer")
    cli::cli_abort("DP must be integer, classification aborted.")

  if(genomic_data(x)$NV %>% class() != "integer")
    cli::cli_abort("NV must be integer, classification aborted.")

  if(any(is.na(clinical_data(x)$purity)) | any(clinical_data(x)$purity < 0 ) | any(!is.numeric(clinical_data(x)$purity)))
    cli::cli_abort("Sample purity must be a non-negative number, classification aborted.")
}

get_sample_priors = function(x, priors, N_mutations, k_max){
  out = lapply(1:nrow(input(x)), function(i){

    gene = input(x)[i,]$gene
    gene_role = input(x)[i,]$gene_role
    tumor_type = input(x)[i,]$tumor_type

    what = priors %>% filter(gene == !!gene, tumor_type == !!tumor_type)

    if(nrow(what) == 0){

      what = priors %>% filter(gene == !!gene, tumor_type == 'PANCA')

      if(nrow(what) == 0){

        if(!is.na(gene_role)){
          what = priors %>% filter(gene == 'other genes', gene_role == !!gene_role, tumor_type == 'PANCA')
        } else {
          what = priors %>% filter(gene == 'other genes', tumor_type == 'PANCA')
        }

      }
    }

    what = what %>% filter(ploidy <= k_max)
    what$p = what$p/sum(what$p)
    what = what %>% arrange(ploidy)
    return(what)
  }) %>% do.call(rbind, .)

  out = lapply(split(out, rep(1:N_mutations, each = k_max*k_max)), function(x){
    c(x$p)
  }) %>% do.call(rbind, .)

}


get_stan_model = function(){
  model_path = system.file("cmdstan", 'model.stan', package = "INCOMMON", mustWork = T)
  model = cmdstanr::cmdstan_model(model_path)
  # tmp = utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
}

get_fit_posterior_per_class = function(fit){
  fit$summary(variables = 'class_probs')
}

get_fit_purity = function(fit){
  fit$summary(variables = 'purity')
}

get_fit_x = function(fit){
  fit$summary(variables = 'x')
}


attach_fit_results = function(x, fit){
  classes = c('m=1','m=k','1<m<k')
  class_probs = get_fit_posterior_per_class(fit)
  outcome = lapply(1:nrow(input(x)), function(i){
    probs = class_probs[grepl(paste0('class_probs\\[',i,','), class_probs$variable),][,'mean']$mean
    names(probs) = classes
    w = probs %>% which.max()
    tibble(
      map_class = classes[w],
      map_posterior = probs[w],
      entropy = -sum(probs * log(probs)),
      probs = list(probs),
      purity_fit = get_fit_purity(fit)[,'median']$median,
      x_fit = get_fit_x(fit)[,'median']$median
    )
  }) %>% do.call(rbind, .)

  dplyr::bind_cols(input(x), outcome)
}


purity = function(x, sample){
  pi = input(x) %>% dplyr::filter(sample == !!sample) %>% dplyr::pull(purity) %>% unique()
  return(pi)
}

samples = function(x){
  input(x) %>% dplyr::pull(sample) %>% unique()
}

# Subset object by sample

subset_sample = function(x, sample){
  stopifnot(inherits(x, 'INCOMMON'))
  samples = unique(x$input$sample)
  stopifnot(sample %in% samples)
  gd = genomic_data(x, PASS = FALSE) %>% dplyr::filter(sample == !!sample)
  cd = clinical_data(x, PASS = FALSE) %>% dplyr::filter(sample == !!sample)
  ip = x$input %>% dplyr::filter(sample == !!sample)
  out = list(genomic_data = gd,
             clinical_data = cd,
             input = ip)
  class(out) = 'INCOMMON'
  if('classification' %in% names(x)) {
    cl = x$classification$fit %>% dplyr::filter(sample == !!sample)
    # pm = x$classification$parameters
    # pr = x$classification$priors


    out$classification$fit = cl
    # out$classification$parameters = pm
    # out$classification$priors = pr
  }
  return(out)
}

#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get classification data for specific selected model.
#' @param x An object of class \code{'INCOMMON'}.
#' @return A table with classified data.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Get classification results
#' classification(MSK_classified)
#' @importFrom dplyr filter mutate rename select %>%
classification = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("classification" %in% names(x))
  stopifnot("fit" %in% names(x$classification))
  x$classification$fit %>%
    # dplyr::full_join(x$input %>%
    #                    dplyr::select(sample, tumor_type, purity) %>% unique(),
    #                  by = 'sample') %>%
    dplyr::select(sample, tumor_type, purity, dplyr::everything())

}
