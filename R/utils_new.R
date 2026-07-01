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

# get_sample_priors = function(x, priors, N_mutations, k_max){
#   out = lapply(1:nrow(input(x)), function(i){
#
#     gene = input(x)[i,]$gene
#     gene_role = input(x)[i,]$gene_role
#     tumor_type = input(x)[i,]$tumor_type
#
#     what = priors %>% filter(gene == !!gene, tumor_type == !!tumor_type)
#
#     if(nrow(what) == 0){
#
#       what = priors %>% filter(gene == !!gene, tumor_type == 'PANCA')
#
#       if(nrow(what) == 0){
#
#         if(!is.na(gene_role)){
#           what = priors %>% filter(gene == 'other genes', gene_role == !!gene_role, tumor_type == 'PANCA')
#         } else {
#           what = priors %>% filter(gene == 'other genes', tumor_type == 'PANCA')
#         }
#
#       }
#     }
#
#     what = what %>% filter(ploidy <= k_max)
#     what$p = what$p/sum(what$p)
#     what = what %>% dplyr::arrange(ploidy)
#     return(what)
#   }) %>% do.call(rbind, .)
#
#   out = lapply(split(out, rep(1:N_mutations, each = k_max*k_max)), function(x){
#     c(x$p)
#   }) %>% do.call(rbind, .)
#
#   return(out)
#
# }

get_sample_priors = function(x, priors, k_max){
  out = lapply(1:nrow(input(x)), function(i){

    gene = input(x)[i,]$gene
    gene_role = input(x)[i,]$gene_role
    tumor_type = input(x)[i,]$tumor_type

    what = priors %>% filter(gene == !!gene, tumor_type == !!tumor_type)

    if(nrow(what) == 0){

      what = priors %>% filter(gene == !!gene, tumor_type == 'PANCA')

      if(nrow(what) == 0){

        if(!is.na(gene_role)){
          what = priors %>%
            dplyr::full_join(INCOMMON::cancer_gene_census) %>%
            dplyr::filter(gene_role == !!gene_role) %>%
            dplyr::group_by(k, m) %>%
            dplyr::reframe(n = mean(n), gene = 'other', gene_role = gene_role, tumor_type = 'PANCA') %>%
            unique() %>%
            dplyr::filter(!is.na(k)) %>%
            dplyr::select(-gene_role)
        } else {
          what = priors %>%
            dplyr::group_by(k, m) %>%
            dplyr::reframe(n = mean(n), gene = 'other', gene_role = NA, tumor_type = 'PANCA') %>%
            unique() %>%
            dplyr::filter(!is.na(k)) %>%
            dplyr::select(-gene_role)
        }

      }
    }

    what = what %>% dplyr::filter(k <= k_max)
    what$N = sum(what$n)
    what$f = what$n/what$N
    what = what %>% dplyr::arrange(k)
    return(what)
  }) %>% do.call(rbind, .)
}

get_stan_input_priors = function(x, priors, N_mutations, k_max){

  out = get_sample_priors(x = x, priors = priors, k_max = k_max)

  out = lapply(split(out, rep(1:N_mutations, each = k_max*(k_max+1)/2)), function(x){
    c(x$n)
  }) %>% do.call(rbind, .)

  return(out)

}


get_stan_model = function(){
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("Package 'cmdstanr' is required to fit INCOMMON models. Install it with:\n  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n  cmdstanr::install_cmdstan()", call. = FALSE)

  model_path = system.file("cmdstan", 'model.stan', package = "INCOMMON", mustWork = T)
  model = cmdstanr::cmdstan_model(model_path)
  model
}

attach_fit_results = function(x, fit, k_max){

  mixing_p_all = fit$summary(variables = 'psi')

  outcome = lapply(1:nrow(input(x)), function(i){

    mixing_p = mixing_p_all[grepl(paste0('psi\\[',i,','), mixing_p_all$variable),][,'mean']$mean

    posterior_table = lapply(1:k_max, function(k){
      lapply(1:k, function(m){
        tibble(k = k,m = m)
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)

    posterior_table$psi = mixing_p
    wc = class_probs %>% which.max()

    k_probs = k_probs[grepl(paste0('posterior_k\\[',i,','), k_probs$variable),][,'mean']$mean
    names(k_probs) = paste0('k=',1:k_max)
    wk = k_probs %>% which.max()

    dplyr::tibble(
      map_k_m = dplyr::arrange(posterior_table, dplyr::desc(psi))[1,],
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

purity_fit = function(x, sample){
  pi = classification(x) %>% dplyr::filter(sample == !!sample) %>% dplyr::pull(purity_fit) %>% unique()
  return(pi)
}

x_fit = function(x, sample){
  x = classification(x) %>% dplyr::filter(sample == !!sample) %>% dplyr::pull(x_fit) %>% unique()
  return(x)
}

samples = function(x){
  input(x) %>% dplyr::pull(sample) %>% unique()
}

# Subset object by sample
#' Subset an INCOMMON object by sample ID.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param sample_list a list of identifiers for the samples to be subsetted
#' @return An object of class `INCOMMON` containing a subset of the original input.
#' @export
#' @examples
#' # First load example data
#' data(MSK_PAAD_output)
#' x = subset_sample(x = MSK_PAAD_output, sample_list = c("P-0000142"))
#' print(x)
#' @importFrom dplyr filter mutate rename select everything %>%
subset_sample = function(x, sample_list){
  stopifnot(inherits(x, 'INCOMMON'))
  samples = unique(x$input$sample)
  stopifnot(length(dplyr::intersect(samples(x), sample_list))>0)
  # gd = genomic_data(x, PASS = FALSE) %>% dplyr::filter(sample %in% sample_list)
  # cd = clinical_data(x, PASS = FALSE) %>% dplyr::filter(sample %in% sample_list)
  x$input = x$input %>% dplyr::filter(sample %in% sample_list)

  return(x)
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

plot_poisson_model = function(x, sample, N_rep, km_rep, km_map, purity_map, eta_map, post_pred_DP, k_max){

  lambda = function(k, x, purity){
    (2*(1-purity)*x+purity*k*x)
  }

  purity_input =  purity(x = x, sample = sample)

  N_stats = get_N_rep_ci(N_rep = N_rep, km_rep = km_rep)

  toplot = x$output %>%
    dplyr::bind_cols(km_map %>% do.call(rbind, .))

  toplot %>%
    dplyr::mutate(test = ifelse((toplot$post_pred_p.value_DP >= .95 | toplot$post_pred_p.value_DP < .05), 'FAIL', 'PASS')) %>%
    ggplot2::ggplot(ggplot2::aes(x = k))+
    ggplot2:: geom_ribbon(data = N_stats, ggplot2::aes(ymin = q5, ymax = q95), linetype=2, alpha=0.5, fill = 'darkgrey')+
    ggplot2::geom_abline(
      data = dplyr::tibble(value = c(purity_map), x_fit = eta_map, purity = c('Input Purity', 'MAP Purity')),
      linetype = 'longdash',
      ggplot2::aes(
        slope = value*eta_map,
        intercept = 2*(1-value)*eta_map), color = 'darkgrey')+
    ggplot2::geom_point(ggplot2::aes(y = DP, fill = test), shape = 21, stroke = 0, size = 3)+
    ggplot2::scale_fill_manual(
      values = c('PASS' = 'forestgreen', 'FAIL' = 'firebrick'))+
    ggrepel::geom_label_repel(ggplot2::aes(x = k, y = DP, label = paste0(gene, ' (',paste(ref, alt, sep = '-'), ')')),
                              force_pull =  -0.05)+
    my_ggplot_theme()+
    ggplot2::scale_x_continuous(expand = c(0.01,0.0))+
    ggplot2::scale_y_continuous(expand = c(0.01,0.0))+
    ggplot2::guides(size = ggplot2::guide_legend(title = 'Posterior Prob'))+
    ggplot2::labs(
      x = 'Total CN (k)',
      y = 'Total reads',
      fill = 'Posterior Predictive test',
      colour = ''
    )+
  ggplot2::xlim(1,max(toplot$k))
}

plot_binomial_model = function(x, n_rep, km_rep, post_pred_NV){

  what = x$output

  what = lapply(1:nrow(what), function(i){
    dplyr::tibble(
      N_rep = n_rep[,i],
      gene = what[i,]$gene,
      id = paste(what[i,]$gene, what[i,]$NV, what[i,]$DP, sep = ':'),
      NV =  what[i,]$NV,
      DP =  what[i,]$DP,
      k = km_rep[[i]]$k,
      m = km_rep[[i]]$m
    )
  }) %>% do.call(rbind, .)

  test = x$output %>%
    dplyr::mutate(id = paste(gene, NV, DP, sep = ':')) %>%
    dplyr::select(id, NV)

  test$p.value = post_pred_NV
  test = test %>% dplyr::mutate(test = ifelse(p.value > .05, 'PASS', 'FAIL'))

  what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(
      ggplot2::aes(
        x = N_rep,
        fill = factor(m)
      ),
      alpha = 0.8,
      binwidth = 1
    ) +
    ggplot2::geom_vline(
      data = test,
      ggplot2::aes(xintercept = NV, group = id, color = test), linetype = 'longdash')+
    ggplot2::scale_color_manual(values = c('PASS' = 'forestgreen', 'FAIL' = 'firebrick'))+
    ggplot2::scale_fill_viridis_d()+
    my_ggplot_theme()+
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.text.y.right= ggplot2::element_text(),
      axis.ticks.y.left = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, colour = 'black'),
      strip.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
    )+
    ggplot2::labs(
      y = '',
      x = 'Reads with the variant',
      fill = 'Multiplicity (m)',
      color = 'Posterior Predictive test'
    )+
    ggplot2::facet_wrap(~id, ncol = 1, strip.position = 'left')
}

bayesian_p_value = function(posterior_rep, prior_rep, prior_type = c("beta", "gamma")) {

  what = rbind(
    dplyr::tibble(
      value = posterior_rep,
      source = 'posterior'
    ),
    dplyr::tibble(
      value = prior_rep,
      source = 'prior'
    )
  )

  prior_type = match.arg(prior_type)

  # Filter data for prior and posterior samples
  prior_data = what %>% dplyr::filter(source == "prior") %>% dplyr::pull(value)
  prior_mean = prior_data %>% mean()
  posterior_data = what %>% dplyr::filter(source == "posterior") %>% dplyr::pull(value)

  # Calculate mean of posterior
  posterior_mean = mean(posterior_data)

  # Calculate the test statistic as absolute difference from prior mean
  if (prior_type == "beta") {
    # Test statistic is the absolute deviation from the prior mean for Beta
    T_prior = abs(prior_data - prior_mean)
    T_posterior = abs(posterior_data - prior_mean)
  } else if (prior_type == "gamma") {
    # Test statistic is the relative deviation from prior mean for Gamma
    T_prior = abs(prior_data - prior_mean) / prior_mean
    T_posterior = abs(posterior_data - prior_mean) / prior_mean
  }

  # Calculate Bayesian p-value as the proportion of times T_posterior >= T_prior
  p_value <- mean(T_posterior >= T_prior)

  # Return results
  return(p_value)
}

posterior_predictive_p_value = function(x, posterior_rep, observed_quantity){
  inp = input(x)
  sapply(1:nrow(inp), function(i){
    mean(abs(posterior_rep[,i] - mean(posterior_rep[,i])) >= abs(inp[[i,observed_quantity]] - mean(posterior_rep[,i])))
  })
}

#' Plot posterior predictive check for eta
#'
#' Visualises the posterior and prior predictive distributions of the
#' expected counts per allele (\eqn{\eta}) and reports the Bayesian p-value
#' from the posterior predictive check.
#'
#' @param posterior_eta_rep Numeric vector of posterior replicated \eqn{\eta} values.
#' @param prior_eta_rep Numeric vector of prior replicated \eqn{\eta} values.
#' @param bayes_p Numeric Bayesian p-value for the posterior predictive check.
#'
#' @return A \code{ggplot2} histogram comparing prior and posterior predictive
#' distributions with vertical median reference lines and the Bayesian p-value.
#'
#' @export
#'
#' @importFrom dplyr tibble mutate
#' @importFrom ggplot2 ggplot geom_histogram aes scale_alpha_manual geom_vline
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme labs geom_text
#' @importFrom ggplot2 element_rect alpha
#' @importFrom scales pretty_breaks
plot_eta_check = function(posterior_eta_rep, prior_eta_rep, bayes_p){

  what = rbind(
    dplyr::tibble(
      value = posterior_eta_rep,
      source = 'posterior'
    ),
    dplyr::tibble(
      value = prior_eta_rep,
      source = 'prior'
    )
  )

  p_value_text = paste0('Bayesian p-value = ', bayes_p %>% format.pval(digits = 2))
  bayes_p_color = ifelse(bayes_p > 0.05, 'forestgreen', 'firebrick')

  p = what %>%
    mutate(source = ifelse(source == 'prior', 'Prior Distribution', 'Posterior Distribution')) %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(ggplot2::aes(x = value, alpha = source), bins = 100, fill = 'darkgrey')+
    ggplot2::scale_alpha_manual(values = c(0.9,0.35))+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(posterior_eta_rep)), color = ggplot2::alpha('darkgrey', alpha = 0.9), linetype = 'longdash')+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(prior_eta_rep)), color = ggplot2::alpha('darkgrey', alpha = 0.35), linetype = 'longdash')+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3))+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
    my_ggplot_theme()+
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = bayes_p_color,
        fill = NA,
        linewidth = 2)
    )+
    ggplot2::labs(
      alpha = '',
      y = '',
      x = bquote('Expected counts per allele'~eta)
    )

  ref_pos = quantile(what$value)['50%']
  text_pos = ifelse(ref_pos > ggplot2::layer_scales(p)$x$range$range[2]/2, 0.0, quantile(what$value)[4])

  p +
    ggplot2::geom_text(
      data = data.frame(), ggplot2::aes(x = max(what$value), y = Inf, label = p_value_text),
      hjust = 1,
      vjust = 2.5,
      color = bayes_p_color, size = 3,
    )
}

#' Plot posterior predictive check for tumour purity
#'
#' Compares prior and posterior predictive distributions of tumour purity
#' (\eqn{\pi}) and reports the Bayesian p-value for the posterior predictive check.
#'
#' @param posterior_purity_rep Numeric vector of posterior replicated purity values.
#' @param prior_purity_rep Numeric vector of prior replicated purity values.
#' @param bayes_p Numeric Bayesian p-value for the posterior predictive check.
#'
#' @return A \code{ggplot2} histogram showing prior and posterior predictive
#' distributions with median reference lines and Bayesian p-value annotation.
#'
#' @export
#'
#' @importFrom dplyr tibble
#' @importFrom ggplot2 ggplot geom_histogram aes scale_alpha_manual geom_vline
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme labs geom_text
#' @importFrom ggplot2 element_rect alpha
#' @importFrom scales pretty_breaks
plot_purity_check = function(posterior_purity_rep, prior_purity_rep, bayes_p){

  what = rbind(
    dplyr::tibble(
      purity = posterior_purity_rep,
      source = 'posterior'
    ),
    dplyr::tibble(
      purity = prior_purity_rep,
      source = 'prior'
    )
  )

  p_value_text = paste0('Bayesian p-value = ', bayes_p %>% format.pval(digits = 2))
  bayes_p_color = ifelse(bayes_p > 0.05, 'forestgreen', 'firebrick')

  p = what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(ggplot2::aes(x = purity, alpha = source), bins = 100, fill = 'darkgrey')+
    ggplot2::scale_alpha_manual(values = c(0.9,0.35))+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(posterior_purity_rep)), color = ggplot2::alpha('darkgrey', alpha = 0.9), linetype = 'longdash')+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(prior_purity_rep)), color = ggplot2::alpha('darkgrey', alpha = 0.35), linetype = 'longdash')+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1))+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
    my_ggplot_theme()+
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = bayes_p_color,
        fill = NA,
        linewidth = 2)
    )+
    ggplot2::labs(
      alpha = '',
      y = '',
      x = bquote('Tumour purity'~pi)
    )

  ref_pos = quantile(what$purity)['50%']
  text_pos = ifelse(ref_pos > ggplot2::layer_scales(p)$x$range$range[2]/2, 0.0, quantile(what$purity)[4])

  p +
    ggplot2::geom_text(
      data = data.frame(), ggplot2::aes(x = max(what$purity), y = Inf, label = p_value_text),
      hjust = 1,
      vjust = 2.5,
      color = bayes_p_color, size = 3,
    )

}

marginal_priors_k = function(x, sample, k_max){
  x = subset_sample(x, sample_list = sample)
  M = nrow(input(x))
  priors = get_sample_priors(x = x, priors = priors_k_m(x), k_max = k_max)
  priors_k = sapply(1:M, function(m){
    sapply(1:k_max, function(k){
      start = ((k-1)*k_max+1)
      end = k_max*k
      priors[m,][start:end] %>% sum()
    })
  }) %>% t() %>% as.data.frame()
  colnames(priors_k) = 1:k_max
  priors_k = dplyr::tibble(priors_k)
  priors_k$gene = classification(x) %>% dplyr::pull(gene)
  priors_k %>%
    tidyr::pivot_longer(cols = colnames(.)[1:k_max], names_to = 'k', values_to = 'p')
}

plot_priors_k = function(x, sample, k_max){
  priors_k = marginal_priors_k(x, sample, k_max)
  priors_k %>%
    dplyr::mutate(k = as.integer(k)) %>%
    ggplot2::ggplot()+
    ggplot2::geom_bar(ggplot2::aes(x = gene, y = p, fill = factor(k)), stat = 'identity')+
    ggplot2::scale_fill_manual(values = ploidy_colors)+
    ggplot2::guides(fill = ggplot2::guide_legend(title = 'k'))+
    my_ggplot_theme()
}

#' Plot joint prior distribution of total copy number and multiplicity
#'
#' Visualises the joint prior distribution of total copy number (\eqn{k})
#' and mutation multiplicity (\eqn{m}) for each mutation in a sample.
#' The most likely configuration for each mutation is highlighted.
#'
#' @param priors_k_m Prior distribution object for joint \eqn{k,m}.
#' @param x An INCOMMON object.
#' @param k_max Integer specifying the maximum total copy number.
#'
#' @return A faceted \code{ggplot2} contour plot showing the joint prior
#' distribution of copy number and multiplicity for each mutation.
#'
#' @export
#'
#' @importFrom dplyr select left_join mutate group_by arrange slice_head
#' @importFrom tidyr complete
#' @importFrom ggplot2 ggplot aes geom_tile geom_contour_filled geom_point
#' @importFrom ggplot2 scale_fill_viridis_c scale_y_continuous guides labs
#' @importFrom ggplot2 theme facet_wrap coord_cartesian element_text
#' @importFrom ggpubr get_legend as_ggplot
plot_prior_k_m = function(priors_k_m, x, k_max){

  what = get_sample_priors(x = x, priors = priors_k_m, k_max = k_max)

  inp = x$output

  what$gene = lapply(1:nrow(inp), function(i){
    rep(inp[i,]$gene, k_max*(k_max+1)/2)
  }) %>% unlist()

  pg = what %>%
    ggplot2::ggplot(ggplot2::aes(x=k,y=m))+
    ggplot2::geom_tile(ggplot2::aes(fill = log10(n)))+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = 'Prior Probability (concentration)', position = 'bottom'))+
    my_ggplot_theme()

  leg = ggpubr::get_legend(pg)
  leg = ggpubr::as_ggplot(leg)

  toplot =  what %>%
    tidyr::complete(
      k = 1:8,
      m = 1:8,
      fill = list(n = 1)  # Or use a small dummy like 1e-40 if log-scaling
    )

 p = toplot %>%
    ggplot2::ggplot(ggplot2::aes(x=k,y=m))+
    ggplot2::geom_contour_filled(ggplot2::aes(z = log10(n)), bins = 50, alpha = 0.85)+
    ggplot2::guides(fill = 'none')+
    ggplot2::coord_cartesian(expand = FALSE, clip = 'off')+
    ggplot2::geom_point(
      data = toplot %>%
        dplyr::left_join(
          x$output %>% dplyr::select(gene, NV)) %>%
        dplyr::mutate(id = paste(gene, NV, sep = ":")) %>%
        dplyr::group_by(id) %>% arrange(desc(f), .by_group = T) %>% slice_head(n=1),
      color = 'firebrick',
      shape = 8)+
    ggplot2::scale_y_continuous(position = 'right') +
    ggplot2::facet_wrap(~id, ncol = 1, strip.position = 'left')+
    my_ggplot_theme()+
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_blank(),panel.grid = ggplot2::element_line(colour = 'black'),
      legend.position = 'bottom',
      axis.text = ggplot2::element_text(),
      axis.title = ggplot2::element_text()
    )+
    ggplot2::labs(x = 'Total CN (k)', y = 'Multiplicity (m)')

 p = p/leg+plot_layout(heights = c(10,1))

 p

}

#' Compute the expectation value of k, m and FAM from the model full posterior distribution.
#' @param x An object of class INCOMMON.
#' @return The same INCOMMON object with annotated expectation values.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' compute_expectations(x = MSK_PAAD_output)
#' @importFrom dplyr mutate group_by summarise full_join select  %>%
#' @importFrom data.table as.data.table
#' @importFrom tidyr unnest
compute_expectations = function(x){

  x$input = x$input %>%
    dplyr::mutate(id = paste(sample, chr, from, to, ref, alt, gene, NV, DP, sep = ':'))

  dt = data.table::as.data.table(x$input)

  dt = tidyr::unnest(dt %>% select(id, z_km), cols = z_km)

  dt = dt %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      FAM = sum((m / k) * z_km),
      exp_m = sum(m * z_km),
      exp_k = sum(k * z_km)
  )

  x$input  = dplyr::full_join(x$input %>% dplyr::select(-z_km), dt)

  return(x)

}

#' Group patients by gene mutant mutant dosage using gene-role specific thresholds.
#' @param x An object of class INCOMMON.
#' @param TSG_low The lower cutoff for mutant dosage classification of tumour suppressor genes.
#' @param TSG_high The upper cutoff for mutant dosage classification of tumour suppressor genes.
#' @param ONC_low The lower cutoff for mutant dosage classification for oncogenes.
#' @param ONC_high The upper cutoff for mutant dosage classification for oncogenes.
#' @return An object of class INCOMMON with new columns reporting mean FAM and class assignment.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' mutant_dosage_classification(MSK_PAAD_output, TSG_low = .25, TSG_high = .75, ONC_low = .33, ONC_high = .66)
#' @importFrom dplyr case_when mutate group_by reframe across everything %>%
mutant_dosage_classification = function(x, TSG_low = .25, TSG_high = .75, ONC_low = .33, ONC_high = .66){
  x = compute_expectations(x)
  x$input = x$input %>%
    dplyr::mutate(class = dplyr::case_when(
      gene_role == 'TSG' & FAM <= TSG_low ~ 'Low Dosage',
      gene_role == 'oncogene' & FAM <= ONC_low ~ 'Low Dosage',

      gene_role == 'TSG' & FAM > TSG_low & FAM < TSG_high ~ 'Balanced Dosage',
      gene_role == 'oncogene' & FAM > ONC_low & FAM < ONC_high ~ 'Balanced Dosage',

      gene_role == 'TSG' & FAM >= TSG_high ~ 'High Dosage',
      gene_role == 'oncogene' & FAM >= ONC_high ~ 'High Dosage'
    ))

  x$input = x$input %>%
    mutate(tmp = paste0(gene, ' with ', class)) %>%
    dplyr::group_by(sample) %>%
    dplyr::reframe(genotype = paste(tmp, collapse = ', '), dplyr::across(dplyr::everything()))

  return(x)
}

#' Get the posterior distribution over (k,m) configurations for a specific mutation
#' @param x An object of class INCOMMON.
#' @param id The id of the mutation (sample:chr:from:to:ref:alt:NV:DP)
#' @return A table with the estimated posterior distribution.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' show_FAM(MSK_PAAD_output, tumor_type = 'PAAD', gene = 'TP53')
#' @importFrom dplyr filter select starts_with %>%
#' @importFrom tidyr unnest
posterior_k_m = function(x, id){
  x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    tidyr::unnest(z_km) %>%
    dplyr::select(id, sample, gene, NV, DP, dplyr::starts_with('purity'), eta_map, m, k, z_km)
}

#' Get the fraction of alleles with the mutation (FAM) values for a gene and cancer type.
#' @param x An object of class INCOMMON.
#' @param tumor_type The tumour type identifier.
#' @param gene The gene name..
#' @return A table with the estimated FAM values
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' show_FAM(MSK_PAAD_output, tumor_type = 'PAAD', gene = 'TP53')
#' @importFrom dplyr filter select %>%
show_FAM = function(x, tumor_type = NULL, gene = NULL){
  if(!("FAM" %in% colnames(x$input))) x = mutant_dosage_classification(x)
  what = x$input
  if(!is.null(tumor_type)){
    what = what %>% dplyr::filter(tumor_type==!!tumor_type)
  }

  if(!is.null(gene)){
    what = what %>% dplyr::filter(gene==!!gene)
  }
  what %>%
    dplyr::select(sample, gene, gene_role, NV, DP, starts_with('purity'), eta_map, FAM, class)
}


draw_samples = function(fit, num_chains, iter_sampling, M, k_max){

  N_rep = fit$draws(variables = 'N_rep') %>% array(dim = c(num_chains * iter_sampling, M))
  n_rep = fit$draws(variables = 'n_rep') %>% array(dim = c(num_chains * iter_sampling, M))

  eta_rep = fit$draws(variables = 'x') %>% array()
  eta_prior_rep = fit$draws(variables = 'x_rep') %>% array()

  purity_rep = fit$draws(variables = 'purity') %>% array()
  purity_prior_rep = fit$draws(variables = 'purity_rep') %>% array()

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% dplyr::arrange(k, m)

  km_idx = fit$draws(variables = 'km_idx') %>% array(dim = c(num_chains * iter_sampling, M))
  km_rep = lapply(1:M, function(i){
    lapply(1:length(km_idx[,i]), function(j){
      k_m_table[km_idx[j,i],]
    }) %>% do.call(rbind, .)
  })

  return(
    list(
      N_rep = N_rep,
      n_rep = n_rep,
      eta_rep = eta_rep,
      eta_prior_rep = eta_prior_rep,
      purity_rep = purity_rep,
      purity_prior_rep = purity_prior_rep,
      km_rep = km_rep
    )
  )
}
