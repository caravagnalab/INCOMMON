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
  model_path = system.file("cmdstan", 'model_v2.stan', package = "INCOMMON", mustWork = T)
  model = cmdstanr::cmdstan_model(model_path)
  # tmp = utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
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
#' subset_sample(x = MSK_PAAD_output, sample_list = c("P-0000142-T01-IM3"))
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

plot_poisson_model = function(x, sample, N_rep, km_rep, km_map, purity_map, x_map, post_pred_DP, k_max){
  lambda = function(k, x, purity){
    (2*(1-purity)*x+purity*k*x)
  }

  inp = input(x)

  purity_input = purity(x = x, sample = sample)

  N_stats = get_N_rep_ci(N_rep = N_rep, km_rep = km_rep)

  toplot = inp %>%
    dplyr::bind_cols(km_map %>% do.call(rbind, .))

  toplot$p.value = post_pred_DP
  toplot = toplot %>% dplyr::mutate(test = ifelse(p.value > .05, 'PASS', 'FAIL'))

  toplot %>%
    ggplot2::ggplot(ggplot2::aes(x = k))+
    ggplot2:: geom_ribbon(data = N_stats, ggplot2::aes(ymin = q5, ymax = q95), linetype=2, alpha=0.5, fill = 'steelblue')+
    ggplot2::geom_abline(
      data = dplyr::tibble(value = c(purity_input, purity_map), x_fit = x_map, purity = c('input', 'fit')),
      linetype = 'longdash',
      ggplot2::aes(
        slope = value*x_map,
        intercept = 2*(1-value)*x_map,
        color = purity))+
    ggplot2::geom_point(ggplot2::aes(y = DP, fill = test), shape = 21, stroke = 0, size = 3)+
    ggplot2::scale_fill_manual(
      values = c('PASS' = 'forestgreen', 'FAIL' = 'firebrick'))+
    ggrepel::geom_label_repel(ggplot2::aes(x = k, y = DP, label = paste(gene, NV, DP, sep = ':')))+
    ggplot2::geom_text(
      data = dplyr::tibble(NULL), label = paste('MAP x =', round(x_map, 2)),
      ggplot2::aes(x = 2, y = N_stats[nrow(N_stats),]$q95))+
    my_ggplot_theme()+
    # ggplot2::ylim(ymin, ymax)+ggplot2::xlim(1, k_max)+
    ggplot2::guides(size = ggplot2::guide_legend(title = 'Posterior Prob'))+
    ggplot2::labs(
      y = 'DP (draws)',
      fill = 'Posterior Predictive test'
      )+
    ggplot2::xlim(1,k_max)
}

plot_binomial_model = function(x, n_rep, km_rep, post_pred_NV){

  what = input(x)
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

  test = input(x) %>%
    dplyr::mutate(id = paste(gene, NV, DP, sep = ':')) %>%
    dplyr::select(id, NV)

  test$p.value = post_pred_NV
  test = test %>% dplyr::mutate(test = ifelse(p.value > .05, 'PASS', 'FAIL'))

  what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(
      ggplot2::aes(
        x = N_rep,
        ),
      fill = 'steelblue',
      alpha = 0.8,
      binwidth = 1
      ) +
    ggplot2::geom_vline(
      # data = what %>%
      #   dplyr::select(id, NV, DP) %>%
      #   unique() %>%
      #   tidyr::pivot_longer(
      #     cols = c('NV', 'DP'),
      #     names_to = 'variable',
      #     values_to = 'value'
      #     ),
      data = test,
      ggplot2::aes(xintercept = NV, color = test, group = id))+
    ggplot2::scale_color_manual(values = c('PASS' = 'forestgreen', 'FAIL' = 'firebrick'))+
    my_ggplot_theme()+
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y.left = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0)
      )+
    ggplot2::labs(
      y = '',
      x = 'NV (draws)',
      fill = 'k:m',
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

plot_x_check = function(posterior_rep, prior_rep, bayes_p){

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

  p_value_text = paste0('Bayesian p-value = ', bayes_p %>% format.pval(digits = 2))
  bayes_p_color = ifelse(bayes_p > 0.05, 'forestgreen', 'firebrick')
  ref_pos = quantile(what$value)['50%'] %>% unname() %>% as.numeric()

  p = what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(ggplot2::aes(x = value, alpha = source), bins = 100, fill = 'steelblue')+
    ggplot2::scale_alpha_manual(values = c(0.8,0.5))+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(posterior_rep)), color = ggplot2::alpha('steelblue', alpha = 0.8), linetype = 'longdash')+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(prior_rep)), color = ggplot2::alpha('steelblue', alpha = 0.5), linetype = 'longdash')+
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=2))+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
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
      x = 'Expected counts per allele'
    )

  text_pos = ifelse(ref_pos > ggplot2::layer_scales(p)$x$range$range[2]/2, 0.0, quantile(what$value)[4])
  p +
    ggplot2::geom_text(
    data = data.frame(), ggplot2::aes(x = text_pos, y = Inf, label = p_value_text),
    hjust = -0.1,
    vjust = 2.5,
    color = bayes_p_color, size = 3,
  )
}

plot_purity_check = function(posterior_rep, prior_rep, bayes_p){

  what = rbind(
    dplyr::tibble(
      purity = posterior_rep,
      source = 'posterior'
    ),
    dplyr::tibble(
      purity = prior_rep,
      source = 'prior'
    )
  )

  p_value_text = paste0('Bayesian p-value = ', bayes_p %>% format.pval(digits = 2))
  bayes_p_color = ifelse(bayes_p > 0.05, 'forestgreen', 'firebrick')
  ref_pos = quantile(what$purity)['50%']
  text_pos = ifelse(ref_pos > .5, 0.0, 0.5)

  what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(ggplot2::aes(x = purity, alpha = source), bins = 100, fill = 'steelblue')+
    ggplot2::scale_alpha_manual(values = c(0.8,0.5))+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(posterior_rep)), color = ggplot2::alpha('steelblue', alpha = 0.8), linetype = 'longdash')+
    ggplot2::geom_vline(ggplot2::aes(xintercept = median(prior_rep)), color = ggplot2::alpha('steelblue', alpha = 0.5), linetype = 'longdash')+
    ggplot2::geom_text(
      data = data.frame(), ggplot2::aes(x = text_pos, y = Inf, label = p_value_text),
      hjust = -.1,
      vjust = 2.5,
      color = bayes_p_color, size = 3,
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=2), limits = c(0,1))+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=2))+
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
      x = 'Sample purity'
    )
}

marginal_priors_k = function(x, sample, k_max){
  x = subset_sample(x, sample_list = sample)
  M = nrow(input(x))
  priors = get_sample_priors(x = x, priors = priors_k_m(x), N_mutations = M, k_max = k_max)
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

plot_prior_k_m = function(priors_k_m, x, k_max){

  what = get_sample_priors(x = x, priors = priors_k_m, k_max = k_max)

  inp = input(x)
  what$gene = lapply(1:nrow(inp), function(i){
    rep(inp[i,]$gene, k_max*(k_max+1)/2)
  }) %>% unlist()

  what %>%
    dplyr::full_join(inp %>% dplyr::select(gene, NV)) %>%
    dplyr::mutate(id = paste0(gene, ':', NV, ' (', tumor_type, ')')) %>%
    ggplot2::ggplot(ggplot2::aes(
    x = factor(k),
    y = factor(m),
    fill = log10(f)
    )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(labels = scales::math_format(10 ^ .x)) +
    ggplot2::scale_x_discrete(breaks = scales::pretty_breaks(n=3))+
    ggplot2::scale_y_discrete(breaks = scales::pretty_breaks(n=3))+
    ggplot2::facet_wrap( ~ id , ncol = 1, strip.position = 'left') +
    my_ggplot_theme() +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, margin = ggplot2::margin())
    )+
    ggplot2::labs(
      x = "Total CN (k)", y = "Multiplicity (m)", fill = "Prior Probability (log10)") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = ggplot2::unit(2.5, 'cm')))
}

plot_km_prior_vs_post = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)
  priors = x$classification$priors_k_m
  k_max = x$classification$parameters$k_max

  priors = get_sample_priors(x = x, priors = priors, k_max = k_max)
  posterior = get_z_km(x = x, sample = sample)

  what = dplyr::bind_cols(
    priors,
    posterior %>% dplyr::select(z_km)
    ) %>%
    tidyr::pivot_longer(cols = c('f', 'z_km'), names_to = 'source', values_to = 'z_km') %>%
    dplyr::mutate(source = ifelse(source == 'f', 'prior', 'posterior'))

  what %>%
    dplyr::filter(source == 'posterior') %>%
    ggplot2::ggplot(ggplot2::aes(
      x = factor(k),
      y = factor(m),
      # fill = round(f, 2)
      fill = log10(z_km)
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_label(ggplot2::aes(label = round(n, 1), color = ''), fill = 'transparent', label.size = 0, size = 3) +
    scale_color_manual(values = c('white'))+
    scale_fill_viridis_c(labels = scales::math_format(10 ^ .x)) +
    my_ggplot_theme() +
    ggh4x::facet_nested_wrap(~tumor_type ~gene)+
    ggplot2::labs(# title = paste0(gene, ' (', tumor_type, ')'),
      x = "Total CN (k)",
      y = "Multiplicity (m)",
      fill = "Posterior Probability (log10)",
      color = 'Dirichlet Prior Concentration'
      ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(barwidth = unit(2.5, 'cm')),
      color = ggplot2::guide_legend(override.aes = list(color = 'white', fill = 'black', size = 5, shape = 21)))

}


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

  return(x)
}

posterior_k_m = function(x, id){
  x = idify(x)
  x$input %>%
    dplyr::filter(id == !!id) %>%
    tidyr::unnest(z_km) %>%
    dplyr::select(id, sample, gene, NV, DP, dplyr::starts_with('purity'), eta_map, m, k, z_km)
}

show_FAM = function(x, tumor_type = NULL, gene = NULL){
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
