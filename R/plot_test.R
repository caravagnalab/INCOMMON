#'Plotting function for class \code{'TAPACLOTH'}.
#' @description
#' Produces a list of plots, one for each mutation in the input, displaying the
#' results of classification.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @import CNAqc
#' @import ggplot2
#' @import ggsci 
#' @import ggridges 
#' @importFrom dplyr filter mutate rename select %>% 
#' @return An object of class \code{'TAPACLOTH'} containing a list of ggplot2
#' plots named `plot_test` inside `classifier`.
#' @export
plot_test = function(x) {
  stopifnot(inherits(x, "TAPACLOTH"))
  x = idify(x)
  plotmodels = lapply(models_avail(x), function(model) {
    plotlist = lapply(get_data(x) %>% pull(id), function(ID) {
      ## Get data for mutation and model used
      mdata = get_classifier(x, model = model) %>% 
        idify() %>% 
        get_data() %>% 
        dplyr::filter(id == ID)
      ## Get sample purity
      purity = get_purity(x)
      ## Format data for plot
      y = lapply(1:(mdata %>% nrow()),
                 function(i) {
                   k = mdata[i,]$karyotype
                   m = mdata[i,]$multiplicity
                   if ((model %>% tolower()) == "binomial") {
                     p = dbinom(
                       x = 1:get_DP(x, ID),
                       size = get_DP(x, ID),
                       prob = m * get_purity(x) / (2 * (1 - get_purity(x)) + get_purity(x) * get_ploidy(k))
                     )
                   }
                   if ((model %>% tolower()) == "beta-binomial") {
                     p = VGAM::dbetabinom(
                       x = 1:get_DP(x, ID),
                       size = get_DP(x, ID),
                       rho = get_rho(x),
                       prob = m * get_purity(x) / (2 * (1 - get_purity(x)) + get_purity(x) * get_ploidy(k))
                     )
                   }
                   tibble(
                     nv = 1:get_DP(x, ID),
                     VAF = nv/unique(mdata$DP),
                     p = p,
                     karyotype = mdata[i,]$karyotype,
                     multiplicity = mdata[i,]$multiplicity,
                     outcome = mdata[i,]$outcome,
                     rank = mdata[i,]$rank,
                     l_a = mdata[i,]$l_a,
                     r_a = mdata[i,]$r_a
                   )
                 }) %>% do.call(rbind, .)
      minrank = y$rank %>% min()
      y = y %>%
        dplyr::mutate(
          class = case_when(
            outcome == "TRUE" & rank == minrank & multiplicity == 1 ~ "BEST1",
            outcome == "TRUE" & rank == minrank & multiplicity == 2 ~ "BEST2",
            outcome == "FALSE" & multiplicity == 1 ~ "FAIL1",
            outcome == "FALSE" & multiplicity == 2 ~ "FAIL2",
            outcome == "TRUE" & multiplicity == 1  & rank != 1 ~ "PASS1",
            outcome == "TRUE" & multiplicity == 2  & rank != 1 ~ "PASS2"
            # outcome == "TRUE" ~ karyotype,
          )
        ) #ifelse(outcome == "FAIL", "FAIL", karyotype))
      ## Build plot
      plt = ggplot2::ggplot() +
        ggridges::geom_density_ridges(
          data = y %>% dplyr::filter(nv >= l_a & nv <= r_a),
          mapping = ggplot2::aes(
            x = nv,
            y = karyotype,
            height = p,
            fill = class,
          ),
          stat = "identity"
        )+scale_x_continuous(sec.axis = sec_axis(trans=~./get_DP(x, ID), name = "VAF"))
      
      if ((model %>% tolower()) == "beta-binomial") {
        model_string = bquote("Test using Beta-Binomial model with " * rho * " = " *
                                .(get_rho(x)))
      } else{
        model_string = bquote("Test using Binomial model")
      }
      
      plt +
        ggplot2::scale_fill_manual(
          values = c(
            "FAIL1" = "#BEBEBE66",
            "FAIL2" = "#BEBEBE66",
            "BEST1" = "firebrick3",
            "BEST2" = "firebrick3",
            "PASS1" =  "#FFA500CC",
            "PASS2" = "#FFA500CC"
          ),
          breaks = c("FAIL", "PASS", "BEST"),
          limits = c(unique(y$class))
          # guide = "none"
        ) + 
        ggplot2::geom_vline(xintercept = get_NV(x, id = ID), linetype = "longdash") +
        CNAqc:::my_ggplot_theme() +
        ggplot2::coord_cartesian(clip = "off", expand = TRUE) +
        ggplot2::labs(
          x = 'NV',
          y = "Pr(X = NV)",
          caption = model_string,
          title = paste0(unique(mdata$ref),">",unique(mdata$alt), "; Gene ", get_gene(x, ID), "(", get_gene_role(x, ID), ")"),
          subtitle = bquote(
            "DP = " * .(unique(mdata$DP)) * '; ' * pi * ' = ' * .(get_purity(x)) * ', Threshold' *
              ' = ' * .(get_threshold(x, model))
          )
        )
    })
    names(plotlist) = get_classifier(x, model) %>% 
      idify() %>% 
      get_data() %>% 
      pull(id) %>% 
      unique()
    return(plotlist)
  })
  names(plotmodels) = models_avail(x)
  for(model in names(plotmodels)) {
    x$classifier[[model]]$plot_test = plotmodels[[model]]
  }
  return(x)
}


