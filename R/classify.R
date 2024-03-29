#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param priors A dplyr::tibble or data frame with columns `gene`, `tumor_type`, `label` and `p` indicating tumor-specific
#' or pan-cancer (PANCA) prior probabilities.
#' @param entropy_cutoff Entropy cut-off for Tier-1 vs Tier-2 assignment.
#' @param rho Over-dispersion parameter.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @return An object of class `INCOMMON` containing the original input plus
#' the classification data and parameters.
#' @export
#' @importFrom dplyr filter mutate rename select %>%

classify = function(x,
                    priors = pcawg_priors,
                    entropy_cutoff = 0.2,
                    rho = 0.01,
                    karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")
)
  {

  stopifnot(inherits(x, "INCOMMON"))
  if(is.null(entropy_cutoff)) entropy_cutoff = 1

  # Output

  output = x

  output$classification = list()
   
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
      w = which(data(x)$id==id)
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

    list(
      fit = dplyr::right_join(input(x) %>% dplyr::select(colnames(genomic_data(x, PASS = TRUE)), id), map, by = 'id'),
      posterior = posterior
    )
  }

    tests = sapply(ids(x), function(id) {
      classify_single_mutation(x = x, id = id)
    })

    output$classification$fit = tests['fit', ] %>% do.call(rbind, .)
    output$classification$posterior = tests['posterior', ]
    output$classification$parameters = dplyr::tibble(entropy_cutoff = entropy_cutoff, rho = rho)

    cli::cli_alert_info('There are: ')
    for (state in c('HMD', 'LOH', 'CNLOH', 'AM', 'Tier-2')) {
      N  = classification(output) %>% dplyr::filter(state == !!state) %>% nrow()
      cli::cli_bullets(c("*" = paste0("N = ", N, ' mutations (', state, ')')))
    }
    
    mean_ent = classification(output) %>% dplyr::pull(entropy) %>% mean()
    min_ent = classification(output) %>% dplyr::pull(entropy) %>% min()
    max_ent = classification(output) %>% dplyr::pull(entropy) %>% max()
    
    cli::cli_alert_info(
      'The mean classification entropy is {.field {round(mean_ent, 2)}} (min: {.field {round(min_ent, 2)}}, max: {.field {round(max_ent, 2)}})'
      )
    
  return(output)
}


