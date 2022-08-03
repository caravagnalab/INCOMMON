#' Compute parameters and distribution for null hypothesis.
#'
#' @param coverage The sequencing coverage at the locus of mutation.
#' @param purity Sample purity.
#' @param rho The overdisperion parameter for "Beta-Binomial" model.
#' @param alpha_level The confidence level of the test.
#' @param model Model used for the test, either "Binomial" or "Beta-Binomial".
#' @return A list.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @import VGAM 
#' @import stringr
#' @examples
#' null_model = test_setup(coverage = 500, purity = 1.0, rho = 0.01, alpha_level = 0.01, model = 'binomial') 
#' print(null_model)
test_setup = function(coverage = 500,
                      purity = 1.0,
                      rho = 0.01,
                      threshold = 0.1,
                      model = 'binomial',
                      karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")
)
{

  # Range of NV values
  nvs = 1:coverage

  # Compute Binomial or Beta-Binomial probability for NV values in range
  log_p = NULL
  # p = purity/2
  # if ((model %>% tolower()) == 'binomial')
  # {
  #   log_p = sapply(nvs, dbinom, size = coverage, prob = p)
  # }
  # else
  # {
  #   log_p = sapply(
  #     nvs,
  #     VGAM::dbetabinom,
  #     size = coverage,
  #     prob = p,
  #     rho = rho
  #   )
  # }
  
  log_p = lapply(karyotypes, function(karyotype) {
    alleles = stringr::str_split(string = karyotype, pattern = ":") %>% unlist() %>% as.integer()
    ploidy = sum(alleles)
    lapply(1:max(alleles), function(mult) {
      if ((model %>% tolower()) == 'binomial') {
        log_p = sapply(nvs,
                       dbinom,
                       size = coverage,
                       prob = mult * purity / (2 * (1 - purity) + purity * ploidy))
      }
      else {
        log_p = sapply(
          nvs,
          VGAM::dbetabinom,
          size = coverage,
          prob = mult * purity / (2 * (1 - purity) + purity * ploidy),
          rho = rho
        )
      }
      
      # Compute P(X > NV) for each NV value
      # p_x = cumsum(log_p)
      # Find l_a such that P(X <= l_a) < alpha
      # l_a = which(p_x < alpha_level, arr.ind = TRUE) %>% max
      # 
      # # Find r_a such that P(X > r_a) < 1 - alpha
      # r_a = which(p_x > 1 - alpha_level, arr.ind = TRUE) %>% min

      # Compute two-tailed p-value for each NV value
      p_x = sapply(nvs, function(nv){log_p[which(log_p <= log_p[nv])] %>% sum()})
      # p_x = 1-p_x
      # p_x = log_p
      
      # Find left extremum NV l_a such that p-value < alpha
      l_a = which(p_x > threshold, arr.ind = TRUE) %>% min
      
      # Find right extremum NV l_a such that p-value < alpha
      r_a = which(p_x > threshold, arr.ind = TRUE) %>% max
      
      # Adjustments for plots when test fails
      if (is.infinite(l_a))
        l_a = 1
      if (is.infinite(r_a))
        r_a = coverage
      
      # Translate NV cutoffs in VAF space
      
      vafs = nvs / coverage
      l_v = vafs[l_a]
      r_v = vafs[r_a]
      
      return(tibble(
        karyotype = karyotype,
        multiplicity = mult,
        inputs = list(tibble(nv = nvs,
                             p = p_x,
                             VAF = vafs)),
          l_a = l_a,
          r_a = r_a,
          l_v = l_v,
          r_v = r_v
      ))
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
  
  
  # cli::cli_alert_info("Computing p-values.")

  return(list(
    model = model,
    rho = rho,
    coverage = coverage,
    purity = purity,
    threshold = threshold,
    test = log_p
  ))
}
