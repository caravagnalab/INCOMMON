# 1. GETTERS
genomic_data = function(x, PASS = TRUE){
  if(PASS) {out = x$genomic_data %>% dplyr::filter(PASS) %>% dplyr::select(-PASS)} else {out = x$genomic_data}
  return(out)
}

clinical_data = function(x, PASS = TRUE){
  if(PASS) {out = x$clinical_data %>% dplyr::filter(PASS) %>% dplyr::select(-PASS)} else {out = x$clinical_data}
  return(out)
}

input = function(x){
  return(x$input)
}

samples = function(x){
  return(x$sample)
}

purity = function(x, id){
  pi = info(x, id = id) %>% dplyr::pull(purity)
  return(pi)
}

tumor_type = function(x, id){
  ttype = info(x, id = id) %>% dplyr::pull(tumor_type)
  return(ttype)
}


#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get model parameters of the performed classification tests.
#' @param x An obj of class \code{'INCOMMON'}.
#' @return A dplyr::tibble containing parameters for all the models used in the classification.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Get classification parameters
#' parameters(MSK_PAAD_output)
#' @importFrom dplyr filter mutate rename select %>%
parameters = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("parameters" %in% names(x))
  x$parameters
}


#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get model priors used for classification.
#' @param x An obj of class \code{'INCOMMON'}.
#' @return A dplyr::tibble containing prior distributions.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Get priors used in classification
#' priors(MSK_PAAD_output)
#' @importFrom dplyr filter mutate rename select %>%
priors_k_m = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("classification" %in% names(x))
  stopifnot("priors_m_k" %in% names(x$classification))
  x$classification$priors_m_k
}

get_categroical_posterior = function(x, id) {

  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("classification" %in% names(x))
  stopifnot("parameters" %in% names(x$classification))

  classification(x) %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(probs) %>% unlist()
}

get_map_posterior = function(x, id) {

  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("classification" %in% names(x))
  stopifnot("parameters" %in% names(x$classification))

  classification(x) %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(map_posterior)
}

idify = function(x){
  x$input = x$input %>%
    dplyr::mutate(id = paste(sample,chr,from,to,ref,alt,NV,DP, sep = ":"))
  return(x)
}

unidify = function(x){
  x$genomic_data = genomic_data(x) %>%
    dplyr::select(-id)
  return(x)
}

ids = function(x){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$input %>% dplyr::pull(id)
}

info = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  out = x$input %>% dplyr::filter(id == !!id)
  if("classification" %in% names(x) & length(x$classification) > 0)
    out = classification(x) %>% dplyr::filter(id == !!id)
  out
}


DP = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(DP)
}

NV = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(NV)
}

entropy = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  posterior(x, id) %>%
    dplyr::filter(NV == NV(x, id)) %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(entropy)
}


VAF = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(VAF)
}

gene = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(gene)
}

get_gene_role = function(x, id){
  if(!("id" %in% colnames(x$input))) x = idify(x)
  x$data %>%
    dplyr::filter(id == !!id) %>%
    dplyr::pull(gene_role)
}


# Prior getter

get_prior = function(x, gene, tumor_type, silent = FALSE){

  if(is.null(x) & !silent) {
    cli::cli_alert("No prior probabilities provided")
    return(1)
  }

  if(!(gene %in% x$gene) & !silent) {
    cli::cli_alert("No prior probability specified for {.field {gene}}")
    return(1)
  }

  if(tumor_type %in% (x %>% dplyr::filter(gene == !!gene) %>% dplyr::pull(tumor_type))) {
    out = x %>% dplyr::filter(gene == !!gene, tumor_type == !!tumor_type)
  } else {
    if(!silent){
      cli::cli_alert("No {.field {tumor_type}}-specific prior probability specified for {.field {gene}}")
      cli::cli_alert("Using a pan-cancer prior")
    }
    out = x %>% dplyr::filter(gene == !!gene, tumor_type == 'PANCA')
    }

  return(out)
}

# Get mutant samples by tumour type and gene

mutant_samples = function(x, tumor_type, gene) {
  stopifnot(inherits(x, 'INCOMMON'))
  input(x) %>%
    dplyr::filter(tumor_type == !!tumor_type,
                  gene == !!gene) %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::reframe(class = unique(class), gene_role) %>%
    dplyr::filter(!grepl('Tier-2', class)) %>%
    dplyr::left_join(clinical_data(x), by = 'sample') %>%
    dplyr::mutate(gene = !!gene)
    # dplyr::left_join(cancer_gene_census, by = 'gene')
}

# Get WT samples by tumour and gene

wt_samples = function(x, tumor_type, gene) {
  stopifnot(inherits(x, 'INCOMMON'))
  input(x) %>%
    dplyr::filter(tumor_type == !!tumor_type) %>%
    dplyr::filter(!grepl(!!gene, genotype)) %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::reframe(class = unique(class)) %>%
    dplyr::mutate(gene = !!gene) %>%
    dplyr::left_join(clinical_data(x), by = 'sample') %>%
    dplyr::mutate(class = paste0(gene, ' WT')) %>%
    dplyr::ungroup()
    # dplyr::left_join(cancer_gene_census, by = 'gene')
}

# 2. CHEKS AND SANITISERS

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

# 3. FORMAT EDITING

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

# 4. ANALYSIS TOOLS


# Get class distribution for an intere cohort

class_frequency = function(x, tumor_type = NULL, gene = NULL, sample_type = NULL){

  what = x$input
  if(!is.null(tumor_type)){
    what = what %>% dplyr::filter(tumor_type==!!tumor_type)
  }

  if(!is.null(gene)){
    what = what %>% dplyr::filter(gene==!!gene)
  }

  if(!is.null(sample_type)){
    what = what %>% dplyr::filter(SAMPLE_TYPE==!!sample_type)
  }

  what %>%
    dplyr::filter(SAMPLE_TYPE!='Unknown') %>%
    dplyr::mutate(SAMPLE_TYPE = ifelse(grepl('Recurrence|Metastasis', SAMPLE_TYPE, ignore.case = T), 'Metastasis', 'Primary')) %>%
    dplyr::select(sample, tumor_type, SAMPLE_TYPE, gene, gene_role, NV, DP, starts_with('purity'), eta_map, FAM, class) %>%
    dplyr::mutate(tot = length(unique(sample))) %>%
    dplyr::group_by(class, tumor_type, SAMPLE_TYPE) %>%
    dplyr::reframe(n = length(unique(sample)), tot = unique(tot)) %>%
    dplyr::group_by(tumor_type, SAMPLE_TYPE) %>%
    dplyr::reframe(N = sum(n), dplyr::across(dplyr::everything())) %>%
    dplyr::mutate(f = n/N)
}


# Prepare input for Kaplan-Meier fit

prepare_km_fit_input = function(x, tumor_type, gene){
  mut = mutant_samples(x = x, tumor_type = tumor_type, gene = gene)
  wt = wt_samples(x = x, tumor_type = tumor_type, gene = gene) %>%
    dplyr::mutate(gene_role = unique(mut$gene_role))

  rbind(
    mut,
    wt
  ) %>%
    dplyr::select('sample', 'tumor_type', 'gene', 'gene_role', dplyr::everything())
}

forest_plot = plot_forest = function(fit, gene, tumor_type, baseline){

  fit$cox_fit[[1]]$call$formula = fit$formula

  cox_fit = fit$cox_fit[[1]]
  s = summary(cox_fit)

  x_limits = c(s$conf.int[,'lower .95'] %>% min(),
               s$conf.int[,'upper .95'] %>% max())

  what = s$conf.int %>% dplyr::as_tibble()
  what$var = rownames(s$conf.int)
  what$p.value = s$coefficients[,ncol(s$coefficients)]

  reference_table =
    lapply(names(cox_fit$xlevels), function(n) {
      dplyr::tibble(
        var = paste0(n, cox_fit$xlevels[n][[1]][1]),
        value = 1,
        low = 1,
        up = 1,
        p.value = NA
      )
    }) %>% do.call(rbind, .)

  leveled = sapply(names(cox_fit$xlevels), function(n){grep(n, what$var)}) %>% unlist() %>% unique()
  not_leveled = what$var[-leveled]
  #
  # reference_table = rbind(
  #   reference_table,
  #   lapply(not_leveled, function(n) {
  #     dplyr::tibble(
  #       var = n,
  #       value = 1,
  #       low = 1,
  #       up = 1,
  #       p.value = NA
  #     )
  #   }) %>% do.call(rbind, .)
  # )

  toplot = dplyr::tibble(
    var = what$var,
    value = what$`exp(coef)`,
    low = what$`lower .95`,
    up = what$`upper .95`,
    p.value = what$p.value
  ) %>%
    rbind(reference_table)

  toplot = toplot %>%
    dplyr::mutate(var = dplyr::case_when(
      grepl('group', var) ~ gsub('group', '', var),
      grepl('sex', var, ignore.case = T) ~ gsub('sex', 'sex: ', var, ignore.case = T),
      grepl('AGE_AT_SEQUENCING', var, ignore.case = T) ~ gsub('AGE_AT_SEQUENCING', 'Age', var, ignore.case = T),
      grepl('TMB_NONSYNONYMOUS', var, ignore.case = T) ~ gsub('TMB_NONSYNONYMOUS', 'TMB', var, ignore.case = T),
      grepl('SAMPLE_TYPE', var, ignore.case = T) ~ gsub('SAMPLE_TYPE', 'Type: ', var, ignore.case = T),
      grepl('^subtype', var, ignore.case = F) ~ gsub('subtype', 'Subtype: ', var, ignore.case = F),
      grepl('patient', var, ignore.case = F) ~ gsub('patient', 'Patient: ', var, ignore.case = F),
      TRUE ~ var
    ))

  what = what %>%
    dplyr::mutate(var = dplyr::case_when(
      grepl('group', var) ~ gsub('group', '', var),
      grepl('sex', var, ignore.case = T) ~ gsub('sex', 'sex: ', var, ignore.case = T),
      grepl('AGE_AT_SEQUENCING', var, ignore.case = T) ~ gsub('AGE_AT_SEQUENCING', 'Age', var, ignore.case = T),
      grepl('TMB_NONSYNONYMOUS', var, ignore.case = T) ~ gsub('TMB_NONSYNONYMOUS', 'TMB', var, ignore.case = T),
      grepl('SAMPLE_TYPE', var, ignore.case = T) ~ gsub('SAMPLE_TYPE', 'Type: ', var, ignore.case = T),
      grepl('^subtype', var, ignore.case = F) ~ gsub('subtype', 'Subtype: ', var, ignore.case = F),
      grepl('patient', var, ignore.case = F) ~ gsub('patient', 'Patient: ', var, ignore.case = F),
      TRUE ~ var
    ))

  levels = c(
    grep('WT', unique(cox_fit$xlevels$group), value = T),
    grep('Mutant', unique(cox_fit$xlevels$group), value = T),
    grep('Low Dosage', unique(cox_fit$xlevels$group), value = T),
    grep('Balanced Dosage', unique(cox_fit$xlevels$group), value = T),
    grep('High Dosage', unique(cox_fit$xlevels$group), value = T)
    # grep('LOH', unique(fit$xlevels$group), value = T),
    # grep('HMD', unique(fit$xlevels$group), value = T)
  )

  # levels = c(levels, setdiff(toplot$var, levels) %>% sort())

  if(any(grepl('sex', what$var))){
    levels = c(levels, 'sex: Female', 'sex: Male')
  }

  if(any(grepl('TMB', what$var))){
    levels = c(levels, "TMB<= 10", "TMB> 10")
  }

  if(any(grepl('Type', what$var))){
    levels = c(levels, "Type: Primary", "Type: Metastasis")
  }

  levels = c(levels, setdiff(toplot$var, levels))

  levels = levels %>% unique()

  toplot = toplot %>%
    dplyr::mutate(var = factor(
      var,
      levels = levels)) %>%
    mutate(p.value = ifelse(is.na(p.value), 1, p.value))

  toplot$var = factor(toplot$var, levels = levels(toplot$var) %>% rev())

  pp = toplot %>%
    dplyr::mutate(num_label = toplot$var %>% seq_along()) %>%
    # dplyr::mutate(num_label = ifelse(var %in% not_leveled & is.na(p.value), NA, num_label)) %>%
    # dplyr::mutate(num_label = ifelse((var %in% not_leveled) & value == 1, NA, num_label)) %>%
    dplyr::mutate(stripe = (num_label%%2==0)) %>%
    dplyr::mutate(id = case_when(
      value == 1 & low == 1 & up  == 1 ~ 'Reference',
      p.value <= .05 ~ 'Significant',
      TRUE ~ 'n.s.'
    )) %>%
    mutate(id = factor(id, levels = c('Reference', 'n.s.', 'Significant'))) %>%
    ggplot2::ggplot(ggplot2::aes(y = var, x = value))+
    # ggplot2::geom_point(ggplot2::aes(color = p.value <= .05))+
    # ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = up, color = p.value <= .05), width = .5)+
    ggplot2::geom_rect(ggplot2::aes(
      ymax = num_label + .5,
      ymin = num_label - .5,
      xmin = -Inf,
      xmax = Inf,
      fill = stripe
    ), alpha = .4)+
    ggplot2::geom_vline(xintercept = 1, linetype = 'longdash', alpha = .5)+
    ggplot2::scale_fill_manual(values = c('gainsboro', 'white'))+
    ggplot2::geom_point(ggplot2::aes(color = id))+
    ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = up, color = id, width = .3))+
    ggplot2::scale_color_manual(values = c('Significant' = 'indianred3', 'n.s.'='black', 'Reference' = 'darkgrey'))+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = x_limits)+
    my_ggplot_theme(cex = .8)+
    ggplot2::ylab('')+
    ggplot2::xlab('Hazard Ratio')+
    ggplot2::guides(fill = 'none')+
    ggplot2::xlim(0, max(toplot$up))+
    ggplot2::labs(color = '')
  pp

}

# Palettes

ploidy_colors = c(
  "steelblue", # Sky Blue
  "forestgreen", # Green
  "#F28E2B",  # Coral
  "#D55E00", # Red
  "navyblue", # Blue
  "deeppink4", # Pink
  "#52BCA3", # Mint Green
  "#999999", # Gray
  "#A3A500", # Olive
  "#FF61C3", # Magenta
  "#E58606", # Dark Orange
  "#5D69B1", # Purple
  "#99C945", # Lime Green
  "#24796C", # Teal
  "#DAA51B", # Mustard
  "#2F8AC4" # Royal Blue
)


# INCOMMON state coloring
scale_color_INCOMMON_class = function(aes = 'fill'){
  colors = c('forestgreen', "pink2", "purple3")
  names(colors) = c("Balanced Dosage","Low Dosage", "High Dosage")
  if(aes == 'fill') {
    ggplot2::scale_fill_manual(values = colors)
  } else {
    ggplot2::scale_color_manual(values = colors)
  }
}

# INCOMMON survival groups

surv_colors = function(baseline = FALSE) {

  if(baseline)
  {out = c('gainsboro', 'steelblue')} else {
      out = c('gainsboro', 'pink2', 'forestgreen','purple3')
  }
  return(out)
}

# INCOMMON higher-level classes

# INCOMMON state coloring
scale_color_INCOMMON_high_level_class = function(aes = 'fill'){
  class_colors = c('Low Dosage' = 'pink2','Balanced Dosage' = 'forestgreen', 'High Dosage' = 'purple3', 'WT' = 'gainsboro', 'ns' = 'gainsboro')
  if(aes == 'fill') {
    ggplot2::scale_fill_manual(values = class_colors)
  } else {
    ggplot2::scale_color_manual(values = class_colors)
  }
}

# Tumour types
scale_color_ttypes = function(aes = 'fill'){
  tumor_colors = c(
    "BLCA" = "#A22E29",
    "BRCA" = "#7D629E",
    "CHOL" = "#B96461",
    "LUAD" = "#CF559A",
    "PAAD" = "#EAB578",
    "CRC" =  "#E4A6A7",
    "UCEC" = "#ABD5D0",
    "PRAD" = "#31482F",
    "MEL" =  "#099668"
  )
  if(aes == 'fill') {
    ggplot2::scale_fill_manual(values = tumor_colors)
  } else {
    ggplot2::scale_color_manual(values = tumor_colors)
  }
}



