#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param priors A tibble or data frame with columns `gene`, `tumor_type`, `label` and `p` indicating tumor-specific
#' or pan-cancer (PANCA) prior probabilities. 
#' @param entropy_cutoff Entropy cut-off for Tier-1 vs Tier-2 assignment.
#' @param rho Over-dispersion parameter.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @return An object of class `INCOMMON` containing the original input plus
#' the classification data and parameters.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' x = classify(
#'     x = x, 
#'     priors = priors,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' print(x)
classify = function(x,
                    priors = NULL,
                    entropy_cutoff = 0.2,
                    rho = 0.01,
                    karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")
                    
)
  {
  stopifnot(inherits(x, "INCOMMON"))
  if(is.null(entropy_cutoff)) entropy_cutoff = 1
  
  # Output
  
  output = x
  
  if (!("fit" %in% names(output)))
    output$fit = list()
  
  cli::cli_h1(
      "INCOMMON inference of copy number and mutation multiplicity for sample {.field {x$sample}}"
    )
    cat("\n")
    
  check_input(x)
    
  cli::cli_alert_info("Performing classification")
    
  x = idify(x)
    
  tests = lapply(ids(x), function(id) {
    # Control for duplicates
    if(info(x, mutation_id = id) %>% nrow() > 1){
      cli_alert_warning(text = "More than one mutation mapped at: {.field {id}}")
      info(x, id)
      
      cli_alert_warning(text = "Keeping first row by default (check your input data)")
      w = which(data(x)$id==id)
      x$data = x$data[-w[2:length(w)],]
    }
    
    # Compute model likelihood, posterior and entropy
    
    compute_posterior(
      NV = NV(x, id),
      DP = DP(x, id),
      gene = gene(x, id),
      priors = priors,
      tumor_type = tumor_type(x),
      purity = purity(x),
      entropy_cutoff = entropy_cutoff,
      rho = rho,
      karyotypes = karyotypes
    )
  })
  
  names(tests) = ids(x)
  
  # Maximum a posteriori classification
  map_estimates = lapply(ids(x), function(id){
    
    map = tests[[id]] %>% 
      dplyr::group_by(NV) %>% 
      dplyr::filter(value == max(value)) %>% 
      dplyr::filter(NV == NV(x, id)) %>% 
      dplyr::ungroup()

    if(nrow(map)>1){
      cli_alert_warning(text =
                          "With purity {.field {purity(x)}} karyotype {.field {map$karyotype}} with multiplicities {.field {map$multiplicity}} have the same likelihood"
                        )
      cli_alert_warning(text =
                          "Simplest case will be selected: {.field {map$karyotype[1]}}"
      )

      map = map[1.,]
    }
    
    map %>% dplyr::select(label, state, value, entropy) %>% dplyr::rename(posterior = value)
}) %>% do.call(rbind, .)
  
  # Add fit object
    output$fit = list(
      classification = dplyr::bind_cols(data(x), map_estimates),
      params = tibble(entropy_cutoff = entropy_cutoff,
                      rho = rho),
      posterior = tests
    )
    
  return(output)
}


