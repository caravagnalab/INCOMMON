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
          what = priors %>% filter(gene == 'other', gene_role == !!gene_role, tumor_type == 'PANCA')
        } else {
          what = priors %>% filter(gene == 'other', tumor_type == 'PANCA') %>%
            group_by(gene, tumor_type, k, m) %>%
            reframe(N = sum(N), n = sum(n), id = paste(NA, tumor_type), gene_role = NA) %>%
            unique()
        }

      }
    }

    what = what %>% filter(k <= k_max)
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
    if(length(intersect(names(x$classification$fit), sample_list))>0){
      cl = x$classification$fit[sample_list]
      pm = x$classification$parameters
      priors_m_k = x$classification$priors_m_k
      priors_x = x$classification$priors_x


      out$classification$fit = cl
      out$classification$parameters = pm
      out$classification$priors_m_k = priors_m_k
      out$classification$priors_x = priors_x

    }
  }
  return(out)
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
  dp = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(DP)
  nv = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(NV)
  x = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(x_fit)
  gene = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(gene)
  tumor_type = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(tumor_type)
  sample = classification(data) %>% dplyr::filter(id == !!id) %>% dplyr::pull(sample)
  purity = purity(x = data, sample = sample)
  k_max = parameters(x = data) %>% dplyr::pull(k_max)

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

plot_poisson_model = function(x, sample, k_max){
  lambda = function(k, x, purity){
    (2*(1-purity)*x+purity*k*x)
  }
  x = subset_sample(x, sample_list = sample)
  purity_fit = classification(x) %>% dplyr::pull(purity_fit) %>% unique()
  purity = purity(x = x, sample = sample)
  x_fit = classification(x) %>% dplyr::pull(x_fit) %>% unique()
  ymin = 2*(1-purity_fit)*x_fit + purity_fit*x_fit
  ymax = 2*(1-purity_fit)*x_fit + purity_fit*x_fit*k_max
  ymin = min(ymin, min(classification(x) %>% dplyr::pull(DP) %>% min))
  ymax = max(ymax, max(classification(x) %>% dplyr::pull(DP) %>% max))

  N_stats = get_N_rep_ci(x = x, sample = sample)

  classification(x) %>%
    dplyr::mutate(k = gsub('k=', '', k)) %>%
    dplyr::mutate(k = as.integer(k)) %>%
    ggplot2::ggplot(ggplot2::aes(x = k))+
    ggplot2:: geom_ribbon(data = N_stats, ggplot2::aes(ymin = q5, ymax = q95), linetype=2, alpha=0.5, fill = 'steelblue')+
    ggplot2::geom_point(ggplot2::aes(y = DP))+
    # ggplot2::geom_point(ggplot2::aes(size = map_k_posterior))+
    # ggplot2::geom_point(ggplot2::aes(x = k, y = expected_dp), color = 'red')+
    ggplot2::geom_abline(
      data = dplyr::tibble(value = c(purity, purity_fit), x_fit = x_fit, purity = c('input', 'fit')),
      linetype = 'longdash',
      ggplot2::aes(
        slope = value*x_fit,
        intercept = 2*(1-value)*x_fit,
        color = purity))+
    ggrepel::geom_label_repel(ggplot2::aes(x = k, y = DP, label = gene))+
    my_ggplot_theme()+
    # ggplot2::ylim(ymin, ymax)+ggplot2::xlim(1, k_max)+
    ggplot2::guides(size = ggplot2::guide_legend(title = 'Posterior Prob'))+
    ggplot2::labs(title = sample,
                  subtitle = paste0('x = ', round(x_fit, 2), '; lambda = ', lambda(k = 1, x = x_fit, purity = purity_fit) %>% round(2)))
}

what

plot_binomial_model = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)
  n_rep = get_n_rep(x = x, sample = sample)
  k_m_rep = get_k_m_rep(x = x, sample = sample)
  niter = x$classification$parameters$stan_stan_iter_sampling * x$classification$parameters$num_chains

  i = 1
  # what = classification(x)[i,]
  what = lapply(1:nrow(classification(x)), function(i){
    tibble(
      N_rep = n_rep[,i],
      gene = classification(x)[i,]$gene,
      id = paste(classification(x)[i,]$gene, classification(x)[i,]$ref, classification(x)[i,]$alt, sep = ':'),
      NV =  classification(x)[i,]$NV,
      DP =  classification(x)[i,]$DP,
      k = k_m_rep[((i-1)*niter+1):(i*niter),]$k,
      m = k_m_rep[((i-1)*niter+1):(i*niter),]$m
    )
  }) %>% do.call(rbind, .)


  what %>%
    ggplot2::ggplot()+
    ggplot2::geom_histogram(ggplot2::aes(x = N_rep, fill = interaction(k, m,sep = ':')), binwidth = 1, alpha = .7)+
    ggplot2::geom_vline(
      data = what %>%
        dplyr::select(id, NV, DP) %>% unique() %>%  tidyr::pivot_longer(cols = c('NV', 'DP'), names_to = 'variable', values_to = 'value'),
      ggplot2::aes(xintercept = value, color = variable), linetype = 'longdash')+
    ggplot2::scale_color_manual(values = c('NV' = 'black', 'DP' = 'firebrick'))+
    # ggplot2::scale_fill_manual()+
    # ggplot2::geom_vline(data = classification(x)[i,], ggplot2::aes(xintercept = DP), linetype = 'longdash')+
    my_ggplot_theme()+
    # ggplot2::xlim(1, classification(x)[i,]$DP)+
    ggplot2::labs(
      y = '',
      fill = 'k:m',
      color = 'Reads'
    )+
    ggplot2::facet_wrap(~id)
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

plot_prior_k_m = function(x, priors, k_max, sample){
  x = subset_sample(x = x, sample_list = sample)
  what = get_sample_priors(x = x, priors = priors, k_max = k_max)
  # what = priors_k_m %>% dplyr::filter(gene == !!gene, tumor_type == !!tumor_type)
  # title = paste0(gene, ' (', tumor_type, ')')
  # if(nrow(what)==0) {
  #   tumor_type = 'P'
  #   what = priors_k_m %>% dplyr::filter(gene == !!gene, tumor_type == !!tumor_type)
  #
  #   gene_role = cancer_gene_census %>% dplyr::filter(gene == !!gene) %>% dplyr::pull(gene_role)
  #   what = priors_k_m %>% dplyr::filter(gene == 'other', gene_role == !!gene_role, tumor_type == !!tumor_type)
  #   title = paste0(gene_role, ' (', tumor_type, ')')
  # }
  ggplot2::ggplot(what, ggplot2::aes(x = factor(k), y = factor(m), fill = round(f,2))) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = round(n,1)), color = "black") +
    ggplot2::scale_fill_gradient(low = "gainsboro", high = "blue", limits = c(min(what$f), max(what$f))) +
    ggplot2::labs(
      # title = paste0(gene, ' (', tumor_type, ')'),
      x = "Total CN (k)",
      y = "Multiplicity (m)",
      fill = "Prior Probability"
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = unit(2.5, 'cm')))+
    my_ggplot_theme()+
    facet_wrap(~id)
}
