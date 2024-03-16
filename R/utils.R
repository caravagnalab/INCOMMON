data = function(x){
  return(x$data)
}

samples = function(x){
  return(x$sample)
}

purity = function(x){
  return(unique(x$purity))
}

tumor_type = function(x){
  return(unique(x$tumor_type))
}


#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get classification data for specific selected model.
#' @param x An object of class \code{'INCOMMON'}.
#' @return A table with classified data.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' x = classify(
#'     x = x, 
#'     priors = pcawg_priors,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' classification(x)
classification = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("fit" %in% names(x))
  stopifnot("classification" %in% names(x$fit))
  x$fit$classification
}

#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get model parameters of the performed classification tests.
#' @param x An obj of class \code{'INCOMMON'}.
#' @return A tibble containing parameters for all the models used in the classification.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' x = classify(
#'     x = x, 
#'     priors = NULL,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' parameters(x)
parameters = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("fit" %in% names(x))
  stopifnot("params" %in% names(x$fit))
  x$fit$params
}


#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get the model posterior distribution of a mutation.
#' @param x An obj of class \code{'INCOMMON'}.
#' @return A table showing posterior distribution and entropy.
#' @export
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' x = classify(
#'     x = x, 
#'     priors = NULL,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' posterior(x, ids(x)[1])
#' 
posterior = function(x, id) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("fit" %in% names(x))
  stopifnot("posterior" %in% names(x$fit))
  stopifnot(id %in% names(x$fit$posterior))
  x$fit$posterior[[id]]
}

idify = function(x){
  x$data = data(x) %>% 
    mutate(id = paste(chr,from,to,ref,alt,sep = ":"))
  return(x)
}

unidify = function(x){
  x$data = data(x) %>% 
    dplyr::select(-id)
  return(x)
}

ids = function(x){
  if(!("id" %in% colnames(data(x)))) x = idify(x)
  data(x) %>% dplyr::pull(id) %>% unique()
}

info = function(x, mutation_id){
  if(!("id" %in% colnames(data(x)))) x = idify(x)
  out = data(x) %>% dplyr::filter(id == mutation_id)
  if("fit" %in% names(x)) out = classification(x) %>% dplyr::filter(id == mutation_id)
  out
}


DP = function(x, id){
  x = idify(x)
  data(x) %>% 
    dplyr::filter(id == !!id) %>% 
    dplyr::pull(DP)
}

NV = function(x, id){
  x = idify(x)
  data(x) %>% 
    dplyr::filter(id == !!id) %>% 
    dplyr::pull(NV)
}


VAF = function(x, mutation_id){
  x = idify(x)
  data(x) %>% 
    dplyr::filter(id == mutation_id) %>% 
    dplyr::pull(VAF)
}

gene = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    dplyr::pull(gene)
}

get_gene_role = function(x, id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == !!id) %>% 
    dplyr::pull(gene_role)
}


# Prior getter

get_prior = function(x, gene, tumor_type){
  
  if(is.null(x)) {
    cli::cli_alert("No prior probabilities provided")
    return(1)
  }
  
  if(!(gene %in% x$gene)) {
    cli::cli_alert("No prior probability specified for {.field {gene}}")
    return(1)
  }
  
  if(tumor_type %in% (x %>% dplyr::filter(gene == !!gene) %>% dplyr::pull(tumor_type))) {
    out = x %>% dplyr::filter(gene == !!gene, tumor_type == !!tumor_type)
  } else {
    cli::cli_alert("No {.field {tumor_type}}-specific prior probability specified for {.field {gene}}")
    cli::cli_alert("Using a pan-cancer prior")
    out = x %>% dplyr::filter(gene == !!gene, tumor_type == 'PANCA')
    } 
  
  return(out)
}


# Input format check

check_input = function(x){
  if(x$sample %>% class() != "character") 
    cli::cli_abort("Unrecogniseable sample id, will not proceed.")
  
  if(x$data$gene %>% class() != "character") 
    cli::cli_abort("Unrecogniseable gene names, will not proceed.")
  
  if(x$data$gene_role %>% class() != "character" | 
     setdiff(x$data$gene_role, c("TSG", "oncogene", NA)) %>% length() != 0) 
    cli::cli_abort("Unrecogniseable gene roles, will not proceed.")
  
  if(x$data$chr %>% class() != "character") 
    cli::cli_abort("Unrecogniseable chromosome names, will not proceed.")
  
  # if(x$data$from %>% class() != "numeric") 
  #   cli::cli_abort("Unrecogniseable mutation start positions \"from\", will not proceed.")
  # 
  # if(x$data$to %>% class() != "numeric") 
  #   cli::cli_abort("Unrecogniseable mutation end positions \"to\", will not proceed.")
  
  if(x$data$ref %>% class() != "character") 
    cli::cli_abort("Unrecogniseable reference alleles, will not proceed.")
  
  if(x$data$alt %>% class() != "character") 
    cli::cli_abort("Unrecogniseable alternative alleles, will not proceed.")
  
  if(x$data$VAF %>% class() != "numeric") 
    cli::cli_abort("Unrecogniseable VAF, will not proceed.")
  
  if(x$data$DP %>% class() != "integer") 
    cli::cli_abort("Unrecogniseable DP, will not proceed.")
  
  if(x$data$NV %>% class() != "integer") 
    cli::cli_abort("Unrecogniseable NV, will not proceed.")
  
  if(is.na(x$purity) | x$purity < 0 | !is.numeric(x$purity))
    cli::cli_abort("Unrecogniseable sample purity, will not proceed.")
}

# Switch to higher-level classification

reduce_classes = function(x) {
  x %>%
    mutate(
      state = case_when(
        label %in% c("2N (Mutated: 1N)") ~ "HMD",
        label %in% c("4N (Mutated: 1N)", "3N (Mutated: 1N)") ~ "LOW MUT/WT",
        label %in% c("1N (Mutated: 1N)") ~ "LOH",
        label %in% c("2N (Mutated: 2N)") ~ "CNLOH",
        label %in% c("3N (Mutated: 2N)", "4N (Mutated: 2N)") ~ "AM"
      )
    )
}

# Extract assignment probability for each higher-level class

add_per_state_probabilities = function(x, NV){
  lapply(1:nrow(x), function(i){
    probs = x[i,]$density[[1]] %>% 
      dplyr::filter(NV == x[i,]$NV) %>% 
      reduce_classes() %>% 
      group_by(state) %>% 
      dplyr::reframe(p_assign = sum(p_assign)) %>% 
      tidyr::pivot_wider(names_from = state, values_from = p_assign)
    tibble(x[i,],
           p_assign_all = list(probs))
  }) %>% do.call(rbind, .)
} 

remove_mutation = function(x, mutation_id){
  x$data = x$data %>% dplyr::filter(id != mutation_id)
  if('fit' %in% names(x)) x$fit = x$fit %>% dplyr::filter(id != mutation_id)
}

# Palettes

ploidy_colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
ploidy_colors = ploidy_colors[c("1:0","1:1","2:1","2:2")]
names(ploidy_colors) = sapply(names(ploidy_colors), function(n){
  strsplit(n,split = ":")[[1]] %>% as.integer() %>% sum()
})
ploidy_colors = c(ploidy_colors, "Tier-2" = 'gray')
