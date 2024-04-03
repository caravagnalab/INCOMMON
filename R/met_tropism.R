#' Logistic regression of metastatic tropism based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' Logistic regression of metastatic tropism based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @param metastic_site The target organ of metastatic diffusion.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit


met_tropism = function(x, gene, tumor_type, metastatic_site){

  stopifnot(inherits(x, 'INCOMMON'))
  stopifnot(metastatic_site %in% clinical_data(x)$METASTATIC_SITE)

  if(!("genotype" %in% (classification(x) %>% names()))) x = genome_interpreter(x)

  data = classification(x) %>%
    dplyr::filter(tumor_type == !!tumor_type,
                  gene == !!gene) %>%
    dplyr::filter(!grepl('Tier-2', class)) %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::reframe(class = unique(class)) %>%
    dplyr::left_join(clinical_data(x), by = 'sample') %>%
    dplyr::filter(SAMPLE_TYPE == 'Metastasis') %>%
    dplyr::mutate(tropic = 1*(METASTATIC_SITE == metastatic_site)) %>%
    dplyr::select(class, METASTATIC_SITE, tropic, SAMPLE_TYPE)

  if(length(unique(data$class))<2 | length(unique(data$tropic))<2) {

    x$metastatic_propensity[[tumor_type]][[gene]] = NULL

    return(x)
  }

  data = data %>%
    dplyr::mutate(class = factor(class)) %>%
    dplyr::mutate(class = stats::relevel(class, ref = grep('without', unique(data$class), value = T)))

  model = stats::glm(data = data, tropic ~ class, family = binomial(link = 'logit'))

  fit = cbind(summary(model)$coefficients, confint(model)) %>% as_tibble()
  fit$var = rownames(confint(model)) %>% gsub('class', '', .)

  fit = dplyr::tibble(
    gene = gene,
    metastatic_site = metastatic_site,
    class = fit[2,]$var,
    OR = fit[2,]$Estimate %>% exp(),
    low = fit[2,]$`2.5 %` %>% exp(),
    up = fit[2,]$`97.5 %` %>% exp(),
    p.value = fit[2,]$`Pr(>|z|)`
  )

  if(!('metastatic_tropism' %in% names(x))) x$metastatic_tropism = list()
  if(!(tumor_type %in% names(x$metastatic_tropism))) x$metastatic_tropism[[tumor_type]] = list()
  if(!(metastatic_site %in% names(x$metastatic_tropism[[tumor_type]]))) x$metastatic_tropism[[tumor_type]][[metastatic_site]] = list()
  if(!(gene %in% names(x$metastatic_tropism[[tumor_type]][[metastatic_site]]))) x$metastatic_tropism[[tumor_type]][[metastatic_site]][[gene]] = list()

  x$metastatic_tropism[[tumor_type]][[metastatic_site]][[gene]] = fit

  print(fit)

  return(x)
}
