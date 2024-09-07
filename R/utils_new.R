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
    what = what %>% dplyr::arrange(ploidy)
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

get_fit_posterior_per_k = function(fit){
  fit$summary(variables = 'posterior_k')
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
  k_probs = get_fit_posterior_per_k(fit)

  outcome = lapply(1:nrow(input(x)), function(i){

    class_probs = class_probs[grepl(paste0('class_probs\\[',i,','), class_probs$variable),][,'mean']$mean
    names(class_probs) = classes
    wc = class_probs %>% which.max()

    k_probs = k_probs[grepl(paste0('posterior_k\\[',i,','), k_probs$variable),][,'mean']$mean
    names(k_probs) = paste0('k=',1:k_max)
    wk = k_probs %>% which.max()

    dplyr::tibble(
      map_class = names(wc),
      map_class_posterior = class_probs[wc],
      entropy = -sum(class_probs * log(class_probs)),
      class_probs = list(class_probs),
      map_k = names(wk),
      map_k_posterior = k_probs[wk],
      entropy_k = -sum(k_probs * log(k_probs)),
      k_probs = list(k_probs),
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

subset_sample = function(x, sample_list){
  stopifnot(inherits(x, 'INCOMMON'))
  samples = unique(x$input$sample)
  stopifnot(length(intersect(samples(x), sample_list))>0)
  gd = genomic_data(x, PASS = FALSE) %>% dplyr::filter(sample %in% sample_list)
  cd = clinical_data(x, PASS = FALSE) %>% dplyr::filter(sample %in% sample_list)
  ip = x$input %>% dplyr::filter(sample %in% sample_list)
  out = list(genomic_data = gd,
             clinical_data = cd,
             input = ip)
  class(out) = 'INCOMMON'
  if('classification' %in% names(x)) {
    if(length(intersect(classification(x)$sample, sample_list))>0){
      cl = x$classification$fit %>% dplyr::filter(sample %in% sample_list)
      pm = x$classification$parameters
      pr = x$classification$priors


      out$classification$fit = cl
      out$classification$parameters = pm
      out$classification$priors = pr
    }
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


my_ggplot_theme = function(cex = 1)
{
  ggplot2::theme_light(base_size = 10 * cex) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.3 * cex, "cm"),
      panel.background = ggplot2::element_rect(fill = "white")
  )
}

get_likelihood = function(dp, nv, m, k, x, purity){

  lambda = function(k, x, purity){
    (2*(1-purity)*x+purity*k*x)
  }

  binom_prob = function(m, k, purity){
    (m*purity)/(2*(1-purity)+k*purity)
  }

  stats::dpois(dp, lambda(k = k, x = x, purity = purity)) * stats::dbinom(x = nv, size = dp, prob = binom_prob(m = m, k = k, purity = purity))
}

compute_likelihood = function(dp, x, purity){
  lapply(1:8, function(k){
    lapply(1:k, function(m){
      nv_x = 1:dp
      value = get_likelihood(dp = dp, nv = nv_x, m = m, k = k, x = x, purity = purity)
      dplyr::tibble(
        k = k,
        m = m,
        nv = nv_x,
        value = value,
        class = dplyr::case_when(
          k == 1 | m == k ~ 'm=k',
          k > 1 & m == 1 ~ 'm=1',
          k > 1 & m > 1 & m < k ~ '1<m<k'
        )
      )
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

#' Plot model likelihood.
#'
#' @param data An object of class \code{'INCOMMON'}.
#' @param id The id of a mutation iin the form sample:chr:from:to:ref:alt:NV:DP.
#' @return A ggplot object.
#' @export
#' @importFrom patchwork wrap_plots plot_annotation
plot_likelihood = function(data, id){
  dp = classification(data) %>% dplyr::filter(id == !!id) %>% pull(DP)
  nv = classification(data) %>% dplyr::filter(id == !!id) %>% pull(NV)
  x = classification(data) %>% dplyr::filter(id == !!id) %>% pull(x_fit)
  gene = classification(data) %>% dplyr::filter(id == !!id) %>% pull(gene)
  tumor_type = classification(data) %>% dplyr::filter(id == !!id) %>% pull(tumor_type)
  sample = classification(data) %>% dplyr::filter(id == !!id) %>% pull(sample)
  purity = purity(x = data, sample = sample)
  k_max = parameters(x = data) %>% pull(k_max)

  likelihood = compute_likelihood(dp = dp, x = x, purity = purity)

  p1 = likelihood %>%
    ggplot2::ggplot(ggplot2::aes(x = nv, y = value, color = factor(k), group = interaction(k,m)))+
    ggplot2::geom_line()+
    ggplot2::geom_point(ggplot2::aes(shape = factor(m)))+
    ggplot2::geom_vline(xintercept = nv, linetype = 'longdash')+
    my_ggplot_theme(cex = .8)+
    ggplot2::guides(color = ggplot2::guide_legend(title = 'Total CN'), shape = ggplot2::guide_legend(title = 'Multiplicity'))

  p2 = likelihood %>%
    dplyr::group_by(class, nv) %>%
    dplyr::reframe(value = sum(value)) %>%
    dplyr::group_by(nv) %>%
    dplyr::reframe(value = value, class) %>%
    ggplot2::ggplot(ggplot2::aes(x = nv, y = value, color = class, group = class))+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept = nv, linetype = 'longdash')+
    my_ggplot_theme(cex = .8)+
    ggplot2::guides(color = ggplot2::guide_legend(title = 'INCOMMON class'))

  patchwork::wrap_plots(p1+p2)+patchwork::plot_annotation(title = id, subtitle = paste0(gene, ' (', tumor_type, ')'))

}
