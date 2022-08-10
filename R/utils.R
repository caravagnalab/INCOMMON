get_data = function(x){
  return(x$data)
}

get_sample = function(x){
  return(x$sample)
}

get_purity = function(x){
  return(unique(x$purity))
}

#' Getter for class \code{'TAPACLOTH'}.
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

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get model parameters of the performed classification tests.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @return A tibble containing parameters for all the models used in the classification.
#' @export
#' 
get_params = function(x, model) {
  stopifnot(inherits(x, "TAPACLOTH"))
    tibble(model = model,
           x$classifier[[model]]$params)
}

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get classification data for specific gene under the selected model.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used in the test from which to get classification data.
#' @param gene_id The name of the gene.
#' @return A tibble with gene-specific classification data.
#' @export
gene_classification = function(x, gene_id, model){
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$classifier[[model]]$data %>% 
    dplyr::filter(gene == gene_id)
  return(y)
}

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get classification data for specific selected model.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used in the test from which to get classification data.
#' @return A tibble with classification data.
#' @export
get_classifier = function(x, model = NULL) {
  stopifnot(inherits(x, "TAPACLOTH"))
  stopifnot("classifier" %in% names(x))
  if (is.null(model)) {
    lapply(names(x$classifier), function(model) {
      x$classifier[[model]]
    }) %>% unlist(recursive = FALSE)
  }
  else{
    model = model %>% tolower()
    x$classifier[[model]]
  }
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

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get coverage for a specific mutation in the sample.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param mutation_id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @return DP of the mutation.
#' @export
get_DP = function(x, id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == !!id) %>% 
    pull(DP)
}

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get number of reads with variant for a specific mutation in the sample.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param mutation_id Coordinates of mutation in the form of a string 
#' containing `chr`,`from`,`to`,`alt`,`ref` coordinates, colon separated.
#' @return NV of the mutation.
#' @export
get_NV = function(x, id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == !!id) %>% 
    pull(NV)
}

#' Getter for class \code{'TAPACLOTH'}.
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

#' Getter for class \code{'TAPACLOTH'}.
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

get_gene_role = function(x, id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == !!id) %>% 
    pull(gene_role)
}

#' Getter for class \code{'TAPACLOTH'}.
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

get_mass = function(x, null_model, id){
  y = null_model$test %>% 
    dplyr::select(karyotype, multiplicity, l_a, r_a)
  
  y$id = id
  
  y$mass = sapply(null_model$test$inputs, function(s) {
    s$p[get_NV(x, id)]
  })
  
  y$gene = get_gene(x, id)
  
  return(y)
}

#' Plot function for class \code{'TAPACLOTH'}.
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
  ids = get_id(x,gene_name)
  lapply(ids, function(i){
    x$classifier[[model]]$plot_test[[i]]
  })
}

#' Plot function for class \code{'TAPACLOTH'}.
#' @description
#' Plot results of the BMix fit used for purity estimation, under the specified model.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimate
#' @return A plot of BMix fit.
#' @export
plot_bmix = function(x, model){
  model = model %>% tolower()
  y = x$purity_estimate[[model]]$plot_bmix
  return(y)
}


#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Extract data and parameters of purity estimation.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param model Model used for purity estimation
#' @return Purity estimate data.
#' @export
get_purity_estimate = function(x, model){
  model = model %>% tolower()
  stopifnot(inherits(x, "TAPACLOTH"))
  y = x$purity_estimate[[model]]
  return(y)
}

#' Getter for class \code{'TAPACLOTH'}.
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

#' Getter for class \code{'TAPACLOTH'}.
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

get_threshold = function(x, model){
  y = get_classifier(x, model)
  y$params$threshold
}

get_rho = function(x){
  y = get_classifier(x, model = "beta-binomial")
  y$params$rho
}

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' List models used for classification tests.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @return Model names.
#' @export
models_avail = function(x){
  stopifnot(inherits(x, "TAPACLOTH"))
  return(names(x$classifier))
}

get_ploidy = function(k){
  stringr::str_split(k, pattern = ":")[[1]] %>% as.integer() %>% sum()}

#' Getter for class \code{'TAPACLOTH'}.
#' @description
#' Get ID of the mutation(s) affecting the specified gene
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param gene_name Name of the selected gene.
#' @return ID( of mutation(s).
#' @export
get_id = function(x, gene_name){
  x = idify(x)
  x$data %>% 
    dplyr::filter(gene == gene_name) %>% 
    pull(id)
}

closer_dist = function(null_model, nv, karyotypes) {
  i = null_model$test %>% 
    rowwise() %>%
    mutate(dist = min(nv - l_a, nv - r_a)) %>% 
    pull(dist) %>% 
    abs() %>% 
    which.min()
  return(i)
}

# Maxima
maximise = function(x)
{
  x %>%
    group_by(NV) %>%
    filter(density == max(density)) %>%
    ungroup()
}
