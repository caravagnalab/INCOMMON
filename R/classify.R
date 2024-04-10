#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param priors A dplyr::tibble or data frame with columns `gene`, `tumor_type`, `label` and `p` indicating tumor-specific
#' or pan-cancer (PANCA) prior probabilities.
#' @param entropy_cutoff Entropy cut-off for Tier-1 vs Tier-2 assignment.
#' @param rho Over-dispersion parameter.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @param parallel Whether to run the classification in parallel (default: FALSE)
#' @param num_cores The number of cores to use for parallel classification.
#' By default, it takes 80% of the available cores.
#' @return An object of class `INCOMMON` containing the original input plus
#' the classification data and parameters.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom parallel mclapply detectCores

classify = function(x,
                    priors = pcawg_priors,
                    entropy_cutoff = 0.2,
                    rho = 0.01,
                    parallel = FALSE,
                    num_cores = NULL,
                    karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")
)
  {

  stopifnot(inherits(x, "INCOMMON"))
  if(is.null(entropy_cutoff)) entropy_cutoff = 1

  cli::cli_h1(
      "INCOMMON inference of copy number and mutation multiplicity for sample {.field {x$sample}}"
    )
    cat("\n")

  check_input(x)

  cli::cli_alert_info("Performing classification")

  x = idify(x)

  classify_single_mutation = function(x, id){

    # Control for duplicates
    if(info(x, id) %>% nrow() > 1){
      cli_alert_warning(text = "More than one mutation mapped at: {.field {id}}")
      info(x, id)

      cli_alert_warning(text = "Keeping first row by default (check your input data)")
      w = which(ids(x)==id)
      x$data = x$data[-w[2:length(w)],]
    }

    # Compute model likelihood, posterior and entropy

    posterior = compute_posterior(
      NV = NV(x, id),
      DP = DP(x, id),
      gene = gene(x, id),
      priors = priors,
      tumor_type = tumor_type(x, id),
      purity = purity(x, id),
      entropy_cutoff = entropy_cutoff,
      rho = rho,
      karyotypes = karyotypes
    )

    # Maximum a posteriori classification

    map = posterior %>%
      dplyr::filter(NV == NV(x, id)) %>%
      dplyr::filter(value == max(value)) %>%
      dplyr::ungroup()

    if (nrow(map) > 1) {
      cli_alert_warning(text =
                          "With purity {.field {purity(x, id)}} karyotype {.field {map$karyotype}} with multiplicities {.field {map$multiplicity}} have the same likelihood")
      cli_alert_warning(text =
                          "Simplest case will be selected: {.field {map$karyotype[1]}}")

      map = map[1., ]
    }

    map = map %>% dplyr::select(label, state, value, entropy) %>% dplyr::rename(posterior = value) %>% dplyr::mutate(id = id)

    fit = dplyr::right_join(input(x) %>% dplyr::select(colnames(genomic_data(x, PASS = TRUE)), id), map, by = 'id')

    return(fit)
  }

  if(parallel){

    if(is.null(num_cores)) num_cores = as.integer(0.8*parallel::detectCores())

    tests = parallel::mclapply(X = ids(x),
                               FUN = classify_single_mutation,
                               x = x,
                               mc.cores = num_cores)
  } else {

    tests = lapply(ids(x), function(id) {
      classify_single_mutation(x = x, id = id)
    })

  }

  # Output

  x$classification = list()

  x$classification$fit = tests %>% do.call(rbind, .)
  x$classification$parameters = dplyr::tibble(
    entropy_cutoff = entropy_cutoff,
    rho = rho,
    karyotypes = list(karyotypes)
  )
  x$classification$priors = priors

    cli::cli_alert_info('There are: ')
    for (state in c('HMD', 'LOH', 'CNLOH', 'AM', 'Tier-2')) {
      N  = classification(x) %>% dplyr::filter(state == !!state) %>% nrow()
      cli::cli_bullets(c("*" = paste0("N = ", N, ' mutations (', state, ')')))
    }

    mean_ent = classification(x) %>% dplyr::pull(entropy) %>% mean()
    min_ent = classification(x) %>% dplyr::pull(entropy) %>% min()
    max_ent = classification(x) %>% dplyr::pull(entropy) %>% max()

    cli::cli_alert_info(
      'The mean classification entropy is {.field {round(mean_ent, 2)}} (min: {.field {round(min_ent, 2)}}, max: {.field {round(max_ent, 2)}})'
      )

  return(x)
}


