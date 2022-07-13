#' Classify mutations per sample using a (Beta-)Binomial model-based test, or
#' per gene using a 3-quantile test.
#'
#' @param data A tibble containing mutations with sample name (sample), gene name
#' (gene), number of reads with variant (nv), coverage (dp), variant allele
#' frequency (VAF), and sample purity (purity) as columns.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model If set to "Binomial" or "Beta-Binomial": classification is run per
#' sample, using a Binomial (no over-dispersion) or Beta-Binomial (over-dispersion
#' included) as expected distribution for the number of reads with variant at
#' fixed coverage and purity. If set to "terzile": classification is run per gene,
#' and is based on a terzile test.
#' @param rho If Beta-Binomial model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#'
#' @return An object of class `TAPACLOTH` that represents the classified input data.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data = list(data = dplyr::tibble(sample = "test", gene = c(paste("test gene ", 1:9), "target gene"), nv = c(seq(10,90,10), 120), dp = c(rep(100, 9), 200), VAF = c(seq(10,90,10), 120)/c(rep(100, 9), 200)), purity = dplyr::tibble(sample = "test", purity = 1))
#' data = run_classifier(x = data, alpha_level = 1e-3, model = "Binomial")
#' print(data)
run_classifier = function(x,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
{
  # Output
  test = list()
  class(test) = "TAPACLOTH"

  stopifnot((model %>% tolower()) %in% c("binomial", "beta-binomial", "terzile"))
  
  if (inherits(x, "TAPACLOTH")) {
    test = x
    x = test$data
    if(!("classifier" %in% names(test))) test$classifier = list()
  }
  else{
    # test$data = x
    test = x
    test$classifier = list()
  }

  if ((model %>% tolower()) %in% c("binomial", "beta-binomial")) {
    
      cli::cli_h1("TAPACLOTH {.field {model}} clonality/Zygosity testing for sample {.field {x$sample}}")
      cat("\n")
      
      cli::cli_alert_info("Computing null model distributions and p-values.")
      
      pvalues = lapply(1:(x$data %>% nrow), function(i) {
        null_model = test_setup(
          coverage = x$data$DP[i],
          purity = x$purity,
          rho = rho,
          alpha_level = alpha_level,
          model = model
        )
        
        # class = run_test(nv = sample_data$nv[i], null_model = null_model)
        # cumprob = null_model$density$p[sample_data$nv[i]]
        pvalues = null_model$test %>% select(karyotype, multiplicity, l_a, r_a)
        pvalues$pvalue = sapply(null_model$test$inputs, function(s) {
          s$p[x$data$NV[i]]
        })
        pvalues$gene = x$data$gene[i]
        pvalues$class = paste0("k=",pvalues$karyotype,",m=",pvalues$multiplicity)
        return(pvalues %>% 
                 dplyr::select(gene, karyotype, multiplicity, class, pvalue,l_a,r_a))
      }) %>% do.call(rbind,.)
      
      # class = sapply(1:(x$data %>% nrow), function(i) {
      #   k = pvalues[[i]]$karyotype[which(pvalues[[i]]$pvalue > alpha_level)]
      #   m = pvalues[[i]]$multiplicity[which(pvalues[[i]]$pvalue > alpha_level)]
      #   paste(k,m,sep = "-",collapse = "|")
      # })
      # 
      # sample_data = sample_data %>%
      #   dplyr::mutate(
      #     p_subclonal = cumprob,
      #     p_loh = 1 - cumprob,
      #     ) %>%
      #   select(-cumprob)
      # 
      # sample_data$p_subclonal = p.adjust(sample_data$p_subclonal, method = "BH")    
      # sample_data$p_loh = p.adjust(sample_data$p_loh, method = "BH")
      # sample_data$class = case_when(
      #   sample_data$p_subclonal <= alpha_level ~ "Subclonal",
      #   sample_data$p_loh <= alpha_level ~ "Clonal LOH",
      #   TRUE ~ "Clonal")
      # 
      # return(sample_data)
    
    # if ((model %>% tolower()) == "binomial") {
    #   test$classifier$binomial = list(
    #     params = tibble(alpha = alpha_level),
    #     data = x %>% select(class, starts_with("p_"))
    #   )
    # }
    # 
    # if ((model %>% tolower()) == "beta-binomial") {
    #   test$classifier$`beta-binomial` = list(
    #     params = tibble(alpha = alpha_level,
    #                   rho = rho),
    #     data = x %>% select(class, starts_with("pvalues"))
    #   )
    # }
      if ((model %>% tolower()) == "binomial") {
        test$classifier$binomial = list(
          params = tibble(alpha = alpha_level),
          data = full_join(test$data, pvalues, by = "gene"))
        # full_join(test$classifier$`beta-binomial`$data, by="gene")
      }
      
      if ((model %>% tolower()) == "beta-binomial") {
        test$classifier$`beta-binomial`= list(
          params = tibble(alpha = alpha_level,
                          rho = rho),
          data = full_join(test$data, pvalues, by = "gene"))
      }
  }

  else{
    # x = lapply(unique(x$data$gene), function(g) {
    #   cli::cli_h1(g)
    # 
    #   gene_data = x$data %>%
    #     dplyr::filter(gene == g)
    # 
    #   terziles = quantile(gene_data$VAF / gene_data$purity, probs = c(0, 1, 0.33)) %>%
    #     round(2)
    # 
    #   gene_data = gene_data %>%
    #     dplyr::mutate(
    #       class = case_when(
    #         VAF / x$puritypurity >= terziles[3] ~ "Clonal LOH",
    #         VAF / purity < terziles[3] ~ "Subclonal/Clonal"
    #       )
    #     )
    # }) %>%
    #   do.call(rbind, .)
    
    x = lapply(unique(x$data$sample), function(s) {
      cli::cli_h1(s)

      sample_data = x$data %>%
        dplyr::filter(sample == s)
      
      sample_purity = dplyr::filter(x$purity, sample == s)$purity

      terziles = quantile(sample_data$VAF / sample_purity, 
                          probs = c(0, 1, 0.33)) %>%
        round(2) %>% sort()

      sample_data = sample_data %>%
        dplyr::mutate(
          class = case_when(
            VAF / sample_purity >= terziles[3] ~ "Clonal LOH",
            VAF / sample_purity < terziles[3] ~ "Subclonal/Clonal"
          )
        )
    }) %>%
      do.call(rbind, .)
    
    test$classifier$terzile = list(
      data = x %>% select(class)
    )
  }

  # test$data = x
  # test$classifier = list(
  #   model = model,
  #   rho = rho,
  #   alpha_level = alpha_level
  # )
  

  return(test)

}
