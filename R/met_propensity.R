#' Logistic regression of metastatic propensity based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' @export
#' @importFrom dplyr filter mutate rename select %>% as_tibble
#' @importFrom survival Surv survfit
#' @importFrom stats glm binomial confint


met_propensity = function(x, gene, tumor_type){
  stopifnot(inherits(x, 'INCOMMON'))
  if(!("genotype" %in% (classification(x) %>% names()))) x = genome_interpreter(x)

  data = classification(x) %>%
    dplyr::filter(tumor_type == !!tumor_type,
                  gene == !!gene) %>%
    dplyr::filter(!grepl('Tier-2', class)) %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::reframe(class = unique(class)) %>%
    dplyr::left_join(clinical_data(x), by = 'sample') %>%
    dplyr::filter(SAMPLE_TYPE == 'Primary') %>%
    dplyr::mutate(metastatic = 1*(MET_COUNT > 0)) %>%
    dplyr::select(class, metastatic)

  if(length(unique(data$class))<2 | length(unique(data$metastatic))<2) {

    x$metastatic_propensity[[tumor_type]][[gene]] = NULL

    return(x)
  }

  data = data %>%
    dplyr::mutate(class = factor(class)) %>%
    dplyr::mutate(class = stats::relevel(class, ref = grep('without', unique(data$class), value = T)))


  model = stats::glm(data = data, metastatic ~ class, family = stats::binomial(link = 'logit'))

  fit = cbind(summary(model)$coefficients, stats::confint(model)) %>% dplyr::as_tibble()
  fit$var = rownames(confint(model)) %>% gsub('class', '', .)

  fit = dplyr::tibble(
    gene = gene,
    class = fit[2,]$var,
    OR = fit[2,]$Estimate %>% exp(),
    low = fit[2,]$`2.5 %` %>% exp(),
    up = fit[2,]$`97.5 %` %>% exp(),
    p.value = fit[2,]$`Pr(>|z|)`
  )

  if(!('metastatic_propensity' %in% names(x))) x$metastatic_propensity = list()
  if(!(tumor_type %in% names(x$metastatic_propensity))) x$metastatic_propensity[[tumor_type]] = list()
  if(!(gene %in% names(x$metastatic_propensity[[tumor_type]]))) x$metastatic_propensity[[tumor_type]][[gene]] = list()

  x$metastatic_propensity[[tumor_type]][[gene]] = fit

  print(fit)

  return(x)
}
