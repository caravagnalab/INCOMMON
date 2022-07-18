get_data = function(x){
  return(x$data)
}

get_sample = function(x){
  return(x$sample)
}

get_purity = function(x){
  return(x$purity)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get significant classification data for the specified model, if already tested.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used in the test from which to get classification data.
#' @return A tibble.
#' @export
get_classes = function(x, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$classifier[[model]]$data %>% 
    dplyr::select(chr, from, to, ref, alt, karyotype, multiplicity, outcome)
  return(y)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get model parameters of the performed classification tests.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used in the test from which to get classification data.
#' @return A tibble containing parameters for all the models used in the classification.
#' @export
#' 
get_params = function(x) {
  stopifnot(inherits(x, "TAPACLOTH"))
  lapply(x$classifier %>% names(), function(model) {
    tibble(model = model,
           x$classifier[[model]]$params)
  }) %>%
    do.call(rbind, .)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get classification data for specific gene under the selected model.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used in the test from which to get classification data.
#' @return A tibble with gene-specific classification data.
#' @export
gene_classification = function(x, gene_id, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$classifier[[model]]$data %>% 
    dplyr::filter(gene == gene_id)
  return(y)
}

get_classifier = function(x, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$classifier[[model]]
  return(y)
}

idify = function(x){
  y = get_data(x)
  y = y %>% 
    mutate(id = paste(chr,from,to,ref,alt,sep = ":"))
  x$data = y
  return(x)
}

unidify = function(x){
  y = get_data(x)
  y = y %>% 
    dplyr::select(-id)
  x$data = y
  return(x)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get coverage for a specific mutation in the sample.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param mutation_id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @return DP of the mutation.
#' @export
get_DP = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(DP)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get number of reads with variant for a specific mutation in the sample.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @param from Start position of the mutation.
#' @return NV of the mutation.
#' @export
get_NV = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(NV)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get variant allele frequency (VAF) for a specific mutation in the sample.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param mutation_id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @return VAF of the mutation.
#' @export
get_VAF = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(VAF)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get ID of the gene affected by the specified mutation.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param mutation_id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @return A tibble with the name of the gene affected by the specified mutation.
#' @export
get_gene = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(gene)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get coordinates of mutation(s) mapped on a specified gene.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param gene_name Name of the gene affected by the mutation.
#' @return A list of mutation coordinates.
#' @export
get_coord = function(x, gene_name){
  x %>% 
    idify() %>% 
    get_data() %>% 
    dplyr::filter(gene == gene_name) %>% 
    pull(id) %>% 
    strsplit(., split = ":")
}

get_pvalues = function(x, null_model, mutation_id){
  y = null_model$test %>% 
    dplyr::select(karyotype, multiplicity, l_a, r_a)
  
  y$pvalue = sapply(null_model$test$inputs, function(s) {
    s$p[get_NV(x, mutation_id)]
  })
  
  y$gene = get_gene(x, mutation_id)
  
  return(y)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Plot results of the classification under the specified model and for the specified gene.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param gene_name Name of gene affected by the mutation.
#' @param model Model used for classification.
#' @return A tibble with classification data for the specified mutation.
#' @export
plot_gene = function(x,model,gene_name){
  stopifnot(inherits(x, "TAPACLOTH"))
  model = model %>% tolower()
  x$classifier[[model]]$plot_test[[get_id(x,gene_name)]]
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Extract data and parameters of purity estimation.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimation
#' @return Purity estimate data.
#' @export
get_purity_estimate = function(x, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$purity_estimate[[model]]
  return(y)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Extract the sample purity as estimated using the specified model
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimation
#' @return Sample purity as obtained by BMix purity estimation procedure.
#' @export
get_purity_bmix = function(x, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = get_purity_estimate(x, model)
  return(y$purity)
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Extract reliability of input sample purity as compared to TAPACLOTH estimate using the specified model.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimation
#' @return Reliability of input purity estimate.
#' @export
get_reliability = function(x, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = get_purity_estimate(x, model)
  return(y$reliability)
}

get_alpha = function(x, model){
  y = get_classifier(x, model)
  y$params$alpha
}

get_rho = function(x){
  y = get_classifier(x, model = "beta-binomial")
  y$params$rho
}

#' #' Getter for class \code{'TAPACLOTH'}.
#' @description
#' List models used for classification tests.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimation
#' @return Model names.
#' @export
models_avail = function(x){
  stopifnot(inherits(x, "TAPACLOTH"))
  return(names(x$classifier))
}

get_ploidy = function(k){
  stringr::str_split(k, pattern = ":")[[1]] %>% as.integer() %>% sum()}
