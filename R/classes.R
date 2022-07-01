classes = function(x, model = "binomial") {
  stopifnot(inherits(x, "TAPACLOTH"))
  if ((model %>% tolower()) %in% c("binomial", "beta-binomial", "terzile"))
  {
    x$data$data %>%
      cbind(x$classifier[[model]]$data) %>%
      as_tibble() %>%
      print()
  }
  else{
    cli::cli_alert("No classification using model {.format {model}}")
  }
}

purity_per_sample = function(x) {
  stopifnot(inherits(x, "TAPACLOTH"))
  
  for (model in x$purity_estimate %>% names()) {
    if ((model %>% tolower()) %in% c("binomial", "beta-binomial"))
    {
      x$data$purity %>%
        dplyr::rename(input_purity = purity) %>%
        full_join(x$purity_estimate[[model]]$reliability, by = "sample") %>%
        full_join(x$purity_estimate[[model]]$purity, by = "sample") %>%
        print()
    }
    
  }
}