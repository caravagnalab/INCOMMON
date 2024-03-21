# 1. GETTERS
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
#' @return A dplyr::tibble containing parameters for all the models used in the classification.
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
#' posterior(x, "chr12:25398285:25398285:C:A")
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

# Get mutant samples by tumour type and gene

mutant_samples = function(x, tumor_type, gene) {
  add_genotypes(x) %>%
    dplyr::filter(tumor_type == !!tumor_type) %>%
    tidyr::separate_rows(genotype, sep = ",\\s*") %>%
    dplyr::filter(grepl(!!gene, genotype, ignore.case = TRUE)) %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::summarise(group = toString(genotype), across(everything())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene = !!gene) %>% 
    dplyr::select(-gene_role) %>% 
    dplyr::left_join(cancer_gene_census, by = 'gene') %>% 
    dplyr::filter(!grepl('Tier-2', group))
}

# Get WT samples by tumour and gene

wt_samples = function(x, tumor_type, gene) {
  add_genotypes(x) %>% 
    dplyr::filter(tumor_type == !!tumor_type) %>%
    dplyr::filter(!grepl(!!gene, genotype)) %>%
    dplyr::mutate(gene = !!gene) %>% 
    dplyr::select(-gene_role) %>% 
    dplyr::left_join(cancer_gene_census, by = 'gene') %>% 
    dplyr::mutate(group = paste0(gene, ' WT'))
}

# 2. CHEKS AND SANITISERS

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

# 3. FORMAT EDITING

# Switch to higher-level classification

reduce_classes = function(x) {
  if(!('state' %in% colnames(x))){
    x = x %>%
      dplyr::mutate(
        state = dplyr::case_when(
          label %in% c("2N (Mutated: 1N)") ~ "HMD",
          label %in% c("4N (Mutated: 1N)", "3N (Mutated: 1N)") ~ "LOW MUT/WT",
          label %in% c("1N (Mutated: 1N)") ~ "LOH",
          label %in% c("2N (Mutated: 2N)") ~ "CNLOH",
          label %in% c("3N (Mutated: 2N)", "4N (Mutated: 2N)") ~ "AM"
        )
      )
  }
 if(!('class' %in% colnames(x)))
   x = x %>%
     dplyr::mutate(
       class = dplyr::case_when(
         gene_role == "TSG" & state %in% c("LOH", "CNLOH") ~ "with LOH",
         gene_role == "oncogene" & state %in% c("AM", "CNLOH") ~ "with AMP",
         gene_role == "TSG" & state == "HMD" ~ "without LOH",
         gene_role == "oncogene" & state == "HMD" ~ "without AMP",
         TRUE ~ "Tier-2"
       )
     )
 return(x)
}

# Add genotypes

add_genotypes = function(x){
  classification_cohort(x) %>%
    reduce_classes() %>%
    dplyr::mutate(group = dplyr::case_when(
      class != 'Tier-2' ~ paste('Mutant', gene, class),
      TRUE ~ paste(class, gene)
    )) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::reframe(genotype = paste(group, collapse = ','), dplyr::across(dplyr::everything()))
}

# 4. ANALYSIS TOOLS

# Get classification data for an intere cohort

classification_cohort = function(x){
  lapply(x, function(x) {
    stopifnot(inherits(x, 'INCOMMON'))
    classification(x) %>% 
      dplyr::mutate(sample = samples(x))
  }) %>% do.call(rbind, .) 
}

# Get class distribution for an intere cohort

class_frequency = function(x, tumor_type, gene){
  frequency_table = classification_cohort(x) %>% 
    dplyr::filter(gene == !!gene, tumor_type == !!tumor_type) %>% 
    dplyr::group_by(state) %>% 
    dplyr::reframe(n = unique(length(sample))) %>% 
    dplyr::mutate(N = sum(n), frequency = n/N)
  frequency_table = cbind(tibble(gene = gene, tumor_type = tumor_type), frequency_table)
  return(frequency_table)
}


# Prepare input for Kaplan-Meier fit

prepare_km_fit_input = function(x, tumor_type, gene){
  rbind(
    mutant_samples(x = x, tumor_type = tumor_type, gene = gene),
    wt_samples(x = x, tumor_type = tumor_type, gene = gene)
  ) %>% 
    dplyr::select('sample', 'tumor_type', 'gene', 'gene_role', dplyr::everything())
}


forest_plot = function(x, tumor_types = FALSE){
  
  if(is.null(x)) return(NULL)
  s = summary(x)
  
  x_limits = c(s$conf.int[,'lower .95'] %>% min(),
               s$conf.int[,'upper .95'] %>% max())
  
  what = s$conf.int %>% as_tibble()
  what$var = rownames(s$conf.int)
  what$p.value = s$coefficients[,ncol(s$coefficients)]
  
  reference_table =
    lapply(names(x$xlevels), function(n) {
      tibble(
        var = paste0(n, x$xlevels[n][[1]][1]),
        value = 1,
        low = 1,
        up = 1,
        p.value = NA
      )
    }) %>% do.call(rbind, .)
  
  toplot = tibble(
    var = what$var,
    value = what$`exp(coef)`,
    low = what$`lower .95`,
    up = what$`upper .95`,
    p.value = what$p.value
  ) %>% 
    rbind(reference_table)
  
  toplot = toplot %>% 
    mutate(var = case_when(
      grepl('group', var) ~ gsub('group', '', var),
      grepl('Sex', var, ignore.case = T) ~ gsub('Sex', 'Sex: ', var, ignore.case = T),
      TRUE ~ var
    )) 
  
  levels = c(x$xlevels$group %>% unique())
  for(c in names(x$xlevels)[-1]){
    levels = c(levels, grep(c, toplot$var, value = T, ignore.case = T) %>% rev())
  }
  
  toplot = toplot %>% 
    dplyr::mutate(var = factor(
      var,
      levels = levels))
  
  toplot$var = factor(toplot$var, levels = levels(toplot$var) %>% rev())
  
  pp = toplot %>% 
    mutate(num_label = toplot$var %>% seq_along()) %>% 
    mutate(stripe = (num_label%%2==0)) %>% 
    ggplot(aes(y = var, x = value))+
    geom_point(aes(color = p.value <= .05))+
    geom_errorbar(aes(xmin = low, xmax = up, color = p.value <= .05), width = .5)+
    geom_rect(aes(
      ymax = num_label + .5,
      ymin = num_label - .5,
      xmin = -Inf,
      xmax = Inf,
      fill = stripe
    ), alpha = .4)+
    geom_vline(xintercept = 1, linetype = 'longdash', alpha = .5)+
    scale_fill_manual(values = c('gainsboro', 'white'))+
    geom_point(aes(color = p.value <= .05))+
    geom_errorbar(aes(xmin = low, xmax = up, color = p.value <= .05), width = .1)+
    scale_color_manual(values = c('indianred3','black') %>% rev())+
    scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = x_limits)+
    CNAqc:::my_ggplot_theme(cex = .8)+
    ylab('')+
    xlab('Hazard Ratio')+
    guides(fill = 'none')
  pp
}

# Palettes

ploidy_colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
ploidy_colors = ploidy_colors[c("1:0","1:1","2:1","2:2")]
names(ploidy_colors) = sapply(names(ploidy_colors), function(n){
  strsplit(n,split = ":")[[1]] %>% as.integer() %>% sum()
})
ploidy_colors = c(ploidy_colors, "Tier-2" = 'gray')

# INCOMMON state coloring 
scale_color_INCOMMON_class = function(aes = 'fill'){
  colors = c("steelblue", "#00468BFF" , "indianred3", "forestgreen", "gainsboro" )
  names(colors) = c("LOH","CNLOH", "AM", "HMD", "Tier-2")
  if(aes == 'fill') {
    scale_fill_manual(values = colors)
  } else {
    scale_color_manual(values = colors)
  }
}

# INCOMMON survival groups

surv_colors = function(gene_role) {
  if(gene_role == 'TSG') 
    out = c('gainsboro', 'forestgreen','goldenrod2')
  else 
    out = c('gainsboro', 'forestgreen','purple3')
  return(out)
}
