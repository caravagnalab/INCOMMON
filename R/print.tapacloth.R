#' #' Plot for class \code{'TAPACLOTH'}.
#'

#' @description
#'
#' If the terzile model is used, the default plot is a histogram of the VAF spectrum is plotted,
#' with expected clonal peak, purity and VAF of mutations affecting the target gene 
#' highlighted. If the Binomial or Beta-Binomial model is used, it includes
#' a plot showing details of the used statistical test in the bottom panel.
#'
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param target_gene The target gene for which classification of mutation calls is wanted.
#' @param sample_name The name of the sample to be analysed.
#' @param... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
plot.TAPACLOTH = function(x, target_gene, sample_name,...){
  
  stopifnot(inherits(x, "TAPACLOTH"))
  
  if(any(class(x$fit) == "bmix"))
  {
    return(x$plot_bmix)
  }
  
  fit_plot = plot_fit(fit = x, target_gene = target_gene, sample_name)
  
  if(x$model == "terzile"){
    return(fit_plot)
  }
  
  else{
    target_data = x$fit %>% 
      dplyr::filter(gene == target_gene)
    
    target_plots = lapply(1:(target_data %>% nrow), function(i){
      
      null_model = test_setup(
        coverage = target_data$dp[i],
        purity = target_data$purity[i],
        rho = x$rho,
        alpha_level = x$alpha_level,
        model = x$model
      )
      
      fit_power = plot_test_power(null_model)+
        ggplot2::geom_vline(xintercept = target_data$nv[i], linetype = 'dashed', size = .5)
      
    }
    )
    
    # Fig assembly
    lp = append(list(fit_plot), target_plots)
    
    figure = ggpubr::ggarrange(
      plotlist = lp,
      nrow = lp %>% length,
      ncol = 1
    )
    
    return(figure)
  }
}

#' #' Print for class \code{'TAPACLOTH'}.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
print.TAPACLOTH = function(x, ...){
  
  stopifnot(inherits(x, "TAPACLOTH"))
  
  if(any(class(x$fit) == "bmix")){
    cli::cli_rule(
      paste(
        crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
        'Purity estimate using BMix: '
      )
    )
    print(x$fit)
    cat("\n")
    cli::cli_rule(
      paste(
        crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
        'Data: '
      )
    )
    print(x$data)
  }
  else{
    if(x$model == "Beta-Binomial"){
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Test using {.field {x$model}} model, with overdispersion parameter {.field {x$rho}} and significance level {.field {x$alpha}}'
        )
      )
    } else{
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Test using {.field {x$model}} model, with significance level {.field {x$alpha}}'
        )
      )
    }
    
    print(x$fit)
  }
  
}
