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
#' data(MSK_classified)
#' # Get classification parameters
#' parameters(MSK_classified)
#' @importFrom dplyr filter mutate rename select %>%
parameters = function(x) {
  stopifnot(inherits(x, "INCOMMON"))
  stopifnot("classification" %in% names(x))
  stopifnot("parameters" %in% names(x$classification))
  x$classification$parameters
}


#' Getter for class \code{'INCOMMON'}.
#' @description
#' Get model priors used for classification.
#' @param x An obj of class \code{'INCOMMON'}.
#' @return A dplyr::tibble containing prior distributions.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Get priors used in classification
#' priors(MSK_classified)
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
  classification(x) %>%
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
  classification(x) %>%
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

# Genome interpreter
#' This function returns, for each mutation, an interpreted state
#' depending on the mutant gene role and its copy-number and multiplicity, among
#' Mutant TSG with or without LOH, and Mutant oncogene with or without amplification.
#' In addition, it assigns an interpreted genome to each sample
#' by integrating the interpreted state of each mutant gene in the sample.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results, as
#' produced by  function `classify`.
#' @return An object or a list of class \code{INCOMMON} with additional columns
#' in `classification`.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Note the outputs to screen
#' genome_interpreter(MSK_classified)
#' @importFrom dplyr filter mutate rename select %>%

# genome_interpreter = function(x){
#   stopifnot(inherits(x, 'INCOMMON'))
#   x$classification$fit = x$classification$fit %>%
#     dplyr::mutate(map_k = as.integer(gsub('k=', '', map_k))) %>%
#     dplyr::mutate(
#       class = dplyr::case_when(
#         gene_role == "TSG" & map_class == 'm=k' ~ paste0("Mutant ", gene, " with LOH"),
#         # gene_role == "oncogene" & (map_class == '1<m<k' | (map_class == 'm=k' & map_k > 1)) ~ paste0("Mutant ", gene, " with AMP"),
#         gene_role == "oncogene" & (map_class == '1<m<k' | map_class == 'm=k') ~ paste0("Mutant ", gene, " with AMP"),
#         gene_role == "TSG" &  map_class == "m=1" ~ paste0("Mutant ", gene, " without LOH"),
#         gene_role == "oncogene" &  map_class == "m=1" ~ paste0("Mutant ", gene, " without AMP"),
#         TRUE ~ paste0(paste0(gene, " Tier-2"))
#       )
#     ) %>% dplyr::group_by(sample) %>%
#     dplyr::reframe(genotype = paste(class, collapse = ','), dplyr::across(dplyr::everything()))
#
#   genotype_table = classification(x) %>%
#     dplyr::group_by(genotype) %>%
#     dplyr::reframe(N = length(unique(sample))) %>%
#     dplyr::arrange(dplyr::desc(N)) %>%
#     dplyr::mutate(frequency = N/sum(N))
#
#   cli::cli_alert_info('There are {.field {nrow(genotype_table)}} different genotypes')
#   cli::cli_alert_info('The most abundant genotypes are:')
#   for(i in 1:3){
#     cli::cli_bullets(c('*' = paste(genotype_table[i,]$genotype,
#                                    paste0("(",genotype_table[i,]$N, " Samples, "),
#                                    paste0("Frequency ",round(genotype_table[i,]$frequency, 2), ")")
#     )
#     )
#     )
#   }
#
#   return(x)
# }
genome_intepreter = function(x, level = 'high'){
  if(level == 'low'){
    classification = lapply(1:nrow(x), function(i){
      x[i,]$z_km[[1]] %>%
        dplyr::mutate(
          map_class = dplyr::case_when(
            k == 1 & k == m ~ 'LOH',
            k > 1 & k == m ~ 'CNLOH',
            k > 2 & m > 1 ~ 'AMP',
            k == 2 & m == 1 ~ 'HMD',
            k > 2 & m == 1 ~ 'Tier-2'
          )) %>%
        dplyr::group_by(map_class) %>%
        dplyr::reframe(z_class = sum(z_km)) %>%
        dplyr::arrange(dplyr::desc(z_class)) %>%
        dplyr::slice_head(n = 1)
    }) %>% do.call(rbind, .)
  } else if(level == 'high'){
    classification = lapply(1:nrow(x), function(i){
      x[i,]$z_km[[1]] %>%
        dplyr::mutate(
          map_class = dplyr::case_when(
            (m / k) <= .5 ~ 'NSA',
            (m / k) > .5 & (m / k) < .8  ~ 'AMP',
            (m / k) >= .8 ~ 'LOH'
          )) %>%
        dplyr::group_by(map_class) %>%
        dplyr::reframe(z_class = sum(z_km)) %>%
        dplyr::arrange(dplyr::desc(z_class)) %>%
        dplyr::slice_head(n = 1)
    }) %>% do.call(rbind, .)
  }
  if('map_class' %in% colnames(x)) x = x %>% select(-c(map_class, z_class))

  x %>%
    dplyr::bind_cols(classification)
}

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

forest_plot = function(x, baseline = FALSE){

  if(is.null(x)) return(NULL)
  s = summary(x)

  x_limits = c(s$conf.int[,'lower .95'] %>% min(),
               s$conf.int[,'upper .95'] %>% max())

  what = s$conf.int %>% dplyr::as_tibble()
  what$var = rownames(s$conf.int)
  what$p.value = s$coefficients[,ncol(s$coefficients)]

  reference_table =
    lapply(names(x$xlevels), function(n) {
      dplyr::tibble(
        var = paste0(n, x$xlevels[n][[1]][1]),
        value = 1,
        low = 1,
        up = 1,
        p.value = NA
      )
    }) %>% do.call(rbind, .)

  leveled = sapply(names(x$xlevels), function(n){grep(n, what$var)}) %>% unlist()
  not_leveled = what$var[-leveled]

  reference_table = rbind(
    reference_table,
    lapply(not_leveled, function(n) {
      dplyr::tibble(
        var = n,
        value = 1,
        low = 1,
        up = 1,
        p.value = NA
      )
    }) %>% do.call(rbind, .)
  )

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
      grepl('Sex', var, ignore.case = T) ~ gsub('Sex', 'Sex: ', var, ignore.case = T),
      TRUE ~ var
    ))

  if(baseline){
    levels = c(
      grep('WT', unique(x$xlevels$group), value = T),
      grep('Mutant', unique(x$xlevels$group), value = T)
    )
  } else {
    levels = c(
      grep('WT', unique(x$xlevels$group), value = T),
      grep('without', unique(x$xlevels$group), value = T),
      grep('with ', unique(x$xlevels$group), value = T)
    )
  }


  for(c in c(names(x$xlevels)[-1], not_leveled)){
    levels = c(levels, grep(c, toplot$var, value = T, fixed = TRUE) %>% rev())
  }

  levels = levels %>% unique()

  toplot = toplot %>%
    dplyr::mutate(var = factor(
      var,
      levels = levels))

  toplot$var = factor(toplot$var, levels = levels(toplot$var) %>% rev())

  pp = toplot %>%
    dplyr::mutate(num_label = toplot$var %>% seq_along()) %>%
    dplyr::mutate(num_label = ifelse(var %in% not_leveled & is.na(p.value), NA, num_label)) %>%
    dplyr::mutate(stripe = (num_label%%2==0)) %>%
    ggplot2::ggplot(ggplot2::aes(y = var, x = value))+
    ggplot2::geom_point(ggplot2::aes(color = p.value <= .05))+
    ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = up, color = p.value <= .05), width = .5)+
    ggplot2::geom_rect(ggplot2::aes(
      ymax = num_label + .5,
      ymin = num_label - .5,
      xmin = -Inf,
      xmax = Inf,
      fill = stripe
    ), alpha = .4)+
    ggplot2::geom_vline(xintercept = 1, linetype = 'longdash', alpha = .5)+
    ggplot2::scale_fill_manual(values = c('gainsboro', 'white'))+
    ggplot2::geom_point(ggplot2::aes(color = p.value <= .05))+
    ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = up, color = p.value <= .05), width = .1)+
    ggplot2::scale_color_manual(values = c('indianred3','black') %>% rev())+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = x_limits)+
    my_ggplot_theme(cex = .8)+
    ggplot2::ylab('')+
    ggplot2::xlab('Hazard Ratio')+
    ggplot2::guides(fill = 'none')
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

surv_colors = function(gene_role, baseline = FALSE) {

  if(baseline)
  {out = c('gainsboro', 'steelblue')} else {
    if(gene_role == 'TSG')
      out = c('gainsboro', 'forestgreen','goldenrod1')
    else
      out = c('gainsboro', 'forestgreen','purple3')
  }
  return(out)
}

# INCOMMON higher-level classes

# INCOMMON state coloring
scale_color_INCOMMON_high_level_class = function(aes = 'fill'){
  class_colors = c('without LOH' = 'forestgreen', 'without AMP' = 'forestgreen', 'with LOH' = 'goldenrod1',
                   'with AMP' = 'purple3', 'Tier-2' = 'gainsboro', 'ns' = 'gainsboro')
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



